#!/bin/env julia

using Helpers, CLI, Printf, DelimitedFiles, Serialization, Variant, Distributions, HypothesisTests, Glob, TSV_old

function somatic(; blacklist="", alignment_path="ar_focus/merged_alignments", protein_altering=true, require_germline=false, ox_filter=true, min_germline_reads=50)
	println("Sample\tChrom\tPosition\tRef allele\tMut allele\tGene\tEffect\tMutant allele%\tWBC allele%\tBG allele%\tGNOMAD\tCOSMIC\tClinVar\tStrand bias\tIGV link")

	# Read mutation blacklist into memory
	blacklist = Set(blacklist == "" ? String[] : [match(r"^\w+\t\w+\t\w+\t\w+", line).match for line in eachline(blacklist)])

	# Read COSMIC database into memory
	cosmic = open(deserialize, "/groups/wyattgrp/users/jbacon/reference/matti_hg38/cosmic_77_hg38.jls")
	gnomad = open(deserialize, "/groups/wyattgrp/users/jbacon/reference/matti_hg38/gnomad_v3.jls")
	clinvar = open(deserialize, "/groups/wyattgrp/users/jbacon/reference/matti_hg38/clinvar-2023-11-12_hg38.jls")

	# Generate a protein effect prediction for all somatic mutations
	variants = String[]
	mut_files = glob("*.mut")
	for line in eachline(pipeline(`cat $mut_files`, `grep Somatic`))
		cols = split(line, '\t')
		push!(variants, join(cols[1:4], '\t'))
	end
	unique!(variants)

	println(stderr, "Predicting effect for $(length(variants)) mutations...")
	effects = predict_effect(variants)
	effects = Dict(variants[k] => effects[k] for k in 1:length(variants))
	variants = []    # Release memory

	for file in glob("*.mut")
		patient = replace(file, ".mut", "")
		println(stderr, patient)
		d = read_tsv(file)
		notes_col = findone(d.headers, "NOTES")
		samples = d.headers[notes_col+1:end-1]
		V = size(d, 1); S = length(samples)

		germline = [r"WBC|gDNA" in s for s in samples]
		without_germline = !any(germline)
		if without_germline
			println(stderr, "Patient $patient has no matched germline samples.")
			if require_germline; continue; end
		end
		for v in 1:V
			if !("Somatic" in d[v, "NOTES"]); continue; end

			chr = String(d[v, 1]); pos = Int32(d[v, 2])

			alt_reads = zeros(Int, S); total_reads = zeros(Int, S)
			starred = falses(S)
			for s in 1:S
				cols = split(d[v, notes_col + s], ':')
				alt_reads[s] = parse(Int, cols[1])
				total_reads[s] = parse(Int, cols[2])
				starred[s] = endswith(d[v, notes_col + s], '*')
			end

			germline_alt = sum(alt_reads[germline])
			germline_total = sum(total_reads[germline])
			if germline_total < min_germline_reads && without_germline == false; continue; end

			mut_signature = join(d[v, 1:4], '\t')
			if mut_signature in blacklist; continue; end
			gene, effect = effects[mut_signature]
			if protein_altering && !is_protein_altering(effect); continue; end

			# Check if the variant is found in GNOMAD
			gnomad_freq = db_value(gnomad, chr, pos, d[v, 3], d[v, 4])
			if gnomad_freq >= 0.005; continue; end

			# Check if the variant is found in COSMIC
			cosmic_score = db_value(cosmic, chr, pos, d[v, 3], d[v, 4])

			clinvar_anno = db_string(clinvar, chr, pos, d[v, 3], d[v, 4])

			for s in findall(starred)
				# Filter based on strand bias
				cols = split(d[v, notes_col + s], ':')
				mut_plus, ref_plus = parse.(Int, split(cols[4], ','))
				if (mut_plus == 0 || mut_plus == 100) && abs(mut_plus - ref_plus) >= 20; continue; end
				if abs(mut_plus - ref_plus) >= 40; continue; end

				# Filter oxidative DNA damage
				if ox_filter && ((d[v, 3] == "C" && d[v, 4] == "A") || (d[v, 3] == "G" && d[v, 4] == "T"))
					min_oxidative_maf = r"OZM-054-900-005.*2018Aug07|GU-19-217.*2020Jan29" in samples[s] ? 0.1 : 0.03
					maf = alt_reads[s] / total_reads[s]
					if maf < min_oxidative_maf; continue; end
				end

				@printf("%s\t", samples[s])
				@printf("%s\t%d\t%s\t%s", d[v, 1], d[v, 2], d[v, 3], d[v, 4])
				@printf("\t%s\t%s", gene, effect)
				@printf("\t%.1f%% (%d / %d)", alt_reads[s] / total_reads[s] * 100, alt_reads[s], total_reads[s])
				if without_germline == false
					@printf("\t%.1f%% (%d / %d)", germline_alt / germline_total * 100, germline_alt, germline_total)
				else
					@printf("\t")
				end
				@printf("\t%.3f%%", d[v, "BACKGROUND%"])
				if gnomad_freq > 0
					@printf("\t%.4f%%", gnomad_freq * 100)
				else
					@printf("\t")
				end
				@printf("\t%s", cosmic_score > 0 ? cosmic_score : "")
				@printf("\t%s", clinvar_anno)

				# Show strand bias as well
				cols = split(d[v, notes_col + s], ':')
				@printf("\t%s%% vs %s%%", mut_plus, ref_plus)

				igv_link = "http://localhost:60151/load?file=http://tambio.ctdna.org/gleas0n/$(alignment_path)/$(samples[s]).bam&genome=hg38&locus=$(d[v, 1]):$(d[v, 2])"
				@printf("\t%s\n", igv_link)
			end
		end
	end
end

mutable struct Mutation
	sample::String
	chrom::String
	position::Int32
	ref_allele::String
	mut_allele::String
	gene::String
	effect::String
	alt_reads::Int32
	total_reads::Int32
	wbc_af::Float32
	bg_af::Float32
	gnomad::Float32
	cosmic::Int32
	clinvar::String
	strand_bias::String
end

function germline_defects(; alignment_path="ar_focus/merged_alignments", min_germline_reads=50, blacklist="")

	# Read mutation blacklist into memory
	blacklist = Set(blacklist == "" ? String[] : [match(r"^\w+\t\w+\t\w+\t\w+", line).match for line in eachline(blacklist)])

	println(stderr, "Reading ClinVar and gnomAD databases into memory...")
	clinvar = open(deserialize, "/groups/wyattgrp/users/jbacon/reference/matti_hg38/clinvar-2023-11-12_hg38.jls")
	gnomad = open(deserialize, "/groups/wyattgrp/users/jbacon/reference/matti_hg38/gnomad_v3.jls")

	germline_variants = Mutation[]
	for file in glob("*.mut")
		patient = replace(file, ".mut", "")
		println(stderr, patient)
		d = read_tsv(file)
		notes_col = findone(d.headers, "NOTES")
		samples = d.headers[notes_col+1:end-1]
		V = size(d, 1); S = length(samples)

		germline = [r"WBC|gDNA|[Bb]enign" in s for s in samples]
		if !any(germline)
			println(stderr, "Patient $patient has no matched germline samples.")
			continue
		end
		for v in 1:V
			chr = String(d[v, 1]); pos = Int32(d[v, 2])

			alt_reads = zeros(Int, S); total_reads = zeros(Int, S)
			for s in 1:S
				cols = split(d[v, notes_col + s], ':')
				alt_reads[s] = parse(Int, cols[1])
				total_reads[s] = parse(Int, cols[2])
			end

			germline_alt = sum(alt_reads[germline])
			germline_total = sum(total_reads[germline])
			if germline_total < 10 || germline_alt / germline_total < 0.15
				continue
			end

			# Check if the variant is found in GNOMAD at > 0.5% frequency
			gnomad_freq = db_value(gnomad, chr, pos, d[v, 3], d[v, 4])
			if gnomad_freq >= 0.005; continue; end

			clinvar_anno = db_string(clinvar, chr, pos, d[v, 3], d[v, 4])

			for s in findall(!germline)
				push!(germline_variants, Mutation(samples[s], chr, pos, d[v, 3], d[v, 4], "", "", alt_reads[s], total_reads[s], germline_alt / germline_total, d[v, "BACKGROUND%"], gnomad_freq, 0, clinvar_anno, ""))
			end
		end
	end

	variants = unique([@sprintf("%s\t%d\t%s\t%s", v.chrom, v.position, v.ref_allele, v.mut_allele) for v in germline_variants])
	effects = Dict(zip(variants, predict_effect(variants)))

	println("Sample\tChrom\tPosition\tRef allele\tMut allele\tGene\tEffect\tSomatic allele%\tGermline allele%\tBG allele%\tGNOMAD\tClinVar\tIGV link")
	for v in germline_variants
		mut_signature = @sprintf("%s\t%d\t%s\t%s", v.chrom, v.position, v.ref_allele, v.mut_allele)
		if mut_signature in blacklist; continue; end

		gene, effect = effects[mut_signature]
		if r"Pathogenic|Likely pathogenic" in v.clinvar || is_truncating(effect)
			# Variant is interesting
		else
			continue
		end

		@printf("%s\t%s\t%s\t%s", v.sample, mut_signature, gene, effect)
		@printf("\t%.1f%% (%d / %d)", v.alt_reads / v.total_reads * 100,
			v.alt_reads, v.total_reads)
		@printf("\t%.1f%%", v.wbc_af * 100)
		@printf("\t%.3f%%", v.bg_af)
		if v.gnomad > 0
			@printf("\t%.4f%%", v.gnomad * 100)
		else
			@printf("\t")
		end
		@printf("\t%s", v.clinvar)

		igv_link = "http://localhost:60151/load?file=http://tambio.ctdna.org/gleas0n/$(alignment_path)/$(v.sample).bam&genome=hg38&locus=$(v.chrom):$(v.position)"
		@printf("\t%s\n", igv_link)
	end
end

function db_value(db, chr, pos, ref, alt)
	chr = drop_chr_prefix(chr)
	anno = get(db.chromosomes, chr, nothing)
	if anno == nothing; return Float32(NaN); end
	matches = searchsorted(anno.position, pos)
	if !isempty(matches)
		change = encode_ref_alt(ref, alt)
		for m in matches
			if anno.change[m] == change
				return Float32(anno.annotation[m])
			end
		end
	end
	return Float32(NaN)
end

function db_string(db, chr, pos, ref, alt)
	chr = drop_chr_prefix(chr)
	anno = get(db.chromosomes, chr, nothing)
	if anno == nothing; return ""; end
	matches = searchsorted(anno.position, pos)
	if !isempty(matches)
		change = encode_ref_alt(ref, alt)
		for m in matches
			if anno.change[m] == change
				return String(anno.annotation[m])
			end
		end
	end
	return ""
end

function predict_effect(variants::Vector{String}; genome="/groups/wyattgrp/users/jbacon/tools/annovar_matti/humandb/hg38")
	genome_version = basename(genome)
	humandb_dir = expanduser(dirname(genome))
	tmp_prefix = @sprintf("annovar-%04X", rand(0:2^16-1))
	reformatted_path = "$tmp_prefix.avinput"

	# Create an input file in ANNOVAR-compatible format. Indels must be changed
	# so they use '-' as the ref or alt allele.
	reformatted_file = open(reformatted_path, "w")	
	for variant in variants
		chr, pos, ref, alt = split(variant, '\t')
		start = parse(Int, pos)
		if length(ref) == 1 && length(alt) > 1
			@assert ref[1] == alt[1]
			ref = "-"; alt = alt[2:end];              # Simple insertion
		elseif length(ref) > 1 && length(alt) == 1
			@assert ref[1] == alt[1]
			ref = ref[2:end]; alt = "-"; start += 1;  # Simple deletion
		elseif length(ref) > 1 && length(alt) > 1
			error("Effect predictions for block substitutions unimplemented!")
		end
		stop = start + length(ref) - 1
		write(reformatted_file, "$chr\t$start\t$stop\t$ref\t$alt\n")
	end
	close(reformatted_file)

	run(pipeline(`/groups/wyattgrp/users/jbacon/tools/annovar_matti/table_annovar.pl $reformatted_path $humandb_dir --buildver $genome_version --remove --outfile $tmp_prefix --operation g --protocol refGene`, stderr=devnull))

	anno_path = "$tmp_prefix.$(genome_version)_multianno.txt"
	anno_file = open(anno_path); readline(anno_file)
	effects = map(eachline(anno_file)) do line
		c = split(line, '\t')
		@assert length(c) >= 7
		gene, effect = translate_annovar_effect(c[7], c[6],
			length(c) >= 9 ? c[9] : "", length(c) >= 10 ? c[10] : "")
		if effect == "Unknown" && !("unknown\tUNKNOWN" in line)
			warn("Cannot parse protein change:\n$line")
		end

		# Special case for TERT promoter mutations
		position = parse(Int, c[2])
		if c[1] == "chr5" && 1_295_069 <= position <= 1_295_568
			gene = "TERT"
			effect = "Promoter mutation"
		end
		return (gene, effect)
	end
	rm(reformatted_path); rm(anno_path)
	return effects
end

mutable struct ChromAnnotations
	position::Vector{Int32}
	change::Vector{Int16}
	#annotation::Vector{Float32}
	annotation::Vector
end

mutable struct Annotations
	prefix::String
	chromosomes::Dict{String, ChromAnnotations}
end

base_num = Dict{String, Int16}("A" => 1, "C" => 2, "G" => 3, "T" => 4)
function encode_ref_alt(ref, alt)
	if length(ref) == 1 && length(alt) == 1
		ref_num = get(base_num, ref, 0)
		alt_num = get(base_num, alt, 0)
		if ref_num == 0 || alt_num == 0; return 0; end
		return ref_num * 16 + alt_num   # Base substitution
	elseif abs(length(ref) - length(alt)) % 3 == 0
		return 1     # Non-frameshift indel
	else
		return 2     # Frameshift indel
	end
end

drop_chr_prefix(chr) = startswith(chr, "chr") ? chr[4:end] : chr

function annotate(db_path::String)
	db = open(deserialize, db_path)
	prefix = db.prefix
	for mut_file in filter(f -> endswith(f, ".mut"), readdir())
		println(mut_file)
		d = readtsv(mut_file, text=true)
		headers = d[1, :]
		notes_col = findone(headers, "NOTES")
		chr = ""; anno = nothing
		for r in 2:size(d, 1)
			if prefix in d[r, notes_col]; continue; end
			if d[r, 1] != chr
				chr = d[r, 1]
				anno = get(db.chromosomes, drop_chr_prefix(chr), nothing)
			end
			if anno == nothing; continue; end

			matches = searchsorted(anno.position, parse(Int32, d[r, 2]))
			if isempty(matches); continue; end

			change = encode_ref_alt(d[r, 3], d[r, 4])
			for m in matches
				if anno.change[m] != change; continue; end
				if anno.annotation[m] isa Float32
					d[r, notes_col] *= @sprintf(" %s:%g.", prefix,
						anno.annotation[m])
				elseif anno.annotation[m] isa String
					d[r, notes_col] *= " $prefix:$(anno.annotation[m])."
				end
			end
		end
		writedlm(mut_file, d)
	end
end

# A mapping from gene names to canonical transcript variant
canonical_transcripts = Set{String}();
genes_with_canonical = Set{String}();

function translate_annovar_effect(gene, func, exonic_func, aa_change)
	if !(r"^((exonic|splicing|intergenic|intronic|upstream|downstream|UTR5|UTR3|ncRNA_exonic|ncRNA_intronic|ncRNA_splicing|ncRNA_UTR5);?)*$" in func)
		error("Unrecognized variant effect '$func'.")
	end

	# Build database of canonical transcripts if not already built
	if isempty(canonical_transcripts)
		canonical = String[]
		genes = String[]
		for r in eachrow(readtsv("/groups/wyattgrp/users/jbacon/reference/matti_hg38/uniprot_canonical_isoforms.tsv")[2:end, :])
			if r[3] == ""; continue; end
			append!(canonical, split(r[3], ','))
			push!(genes, r[1])
		end
		global canonical_transcripts = Set(canonical)
		global genes_with_canonical = Set(genes)
	end

	effects = String[]
	if "ncRNA_exonic" in func
		push!(effects, "ncRNA exon")  # Check this before "exonic" below
	elseif "exonic" in func
		canonical = falses(0)
		for effect in split(aa_change, r"[;,] ?")
			m = match(r"^.+:([NX]M_\d+):.*:p.(.+)", effect)
			if m == nothing; return (gene, "Unknown"); end
			transcript = m.captures[1]
			change = m.captures[2]

			m = match(r"^([A-Z*])\d+(.+)|(\d+_\d+)(del)", change)
			if m == nothing; return (gene, "Unknown"); end
			
			if m.captures[1] == m.captures[2]
				push!(effects, "Synonymous p.$(change)")
			elseif "fs" in m.captures[2]
				push!(effects, "Frameshift p.$(change)")
			elseif "del" in m.captures[2]
				push!(effects, "Non-frameshift deletion p.$(change)")
			elseif "ins" in m.captures[2]
				push!(effects, "Non-frameshift insertion p.$(change)")
			elseif m.captures[2] == "*"
				push!(effects, "Stopgain p.$(change[1:end-1])X")
			elseif m.captures[2] == "X"
				push!(effects, "Stopgain p.$(change)")
			elseif m.captures[2] == "?"
				push!(effects, "Startloss p.$(change)")
			elseif m.captures[1] == "X" || m.captures[1] == "*"
				push!(effects, "Stoploss p.$(change)")
			elseif length(m.captures[1]) == length(m.captures[2]) == 1
				push!(effects, "Missense p.$(change)")
			else
				return (gene, "Unknown")
			end

			push!(canonical, transcript in canonical_transcripts)
		end

		# Only keep canonical isoform effects, if there are any
		if any(canonical)
			effects = effects[canonical]
		elseif gene in genes_with_canonical
			push!(effects, "Only affects a non-canonical isoform")
		end
		effects = unique(effects)

	elseif "UTR3" in func; push!(effects, "3'-UTR");
	elseif "UTR5" in func; push!(effects, "5'-UTR");
	elseif "upstream" in func; push!(effects, "Upstream");
	elseif "downstream" in func; push!(effects, "Downstream");
	elseif "intronic" in func; push!(effects, "Intronic");
	elseif "intergenic" in func || func == ""
		push!(effects, "Intergenic")
		gene = ""    # TODO: This line needed?
	end

	if "splicing" in func; pushfirst!(effects, "Splice site"); end

	if isempty(effects); return (gene, "Unknown"); end


	return (gene, join(effects, ". "))
end

subcommands(somatic, germline_defects, predict_effect, annotate)

