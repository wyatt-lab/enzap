
module Variant

using Helpers, DelimitedFiles, Distributions, Printf

export VCF, MinimalVCF
export read_vcf, read_minimal_vcf, read_vcf_sample_names
export is_protein_altering, is_missense, is_truncating, is_indel, is_coding_region
export truncal_maf_to_cancer_frac_assuming_loh
export gene_mutated, mutation_table

mutable struct VCF
	sample::Vector{String}
	chromosome::Vector{String}
	position::Vector{Int32}
	ref_allele::Vector{String}
	alt_allele::Vector{String}
	gene::Vector{String}
	effect::Vector{String}
	notes::Vector{String}
	alt::Matrix{Int32}
	total::Matrix{Int32}
	star::BitMatrix
	mapq::Matrix{Int32}
	sidedness::Matrix{Int32}
end
VCF() = VCF([], [], [], [], [], [], [], [], zeros(0, 0), zeros(0, 0),
	falses(0, 0), zeros(0, 0), zeros(0, 0))
Base.getindex(vcf::VCF, rows, cols) = VCF(vcf.sample[cols],
	vcf.chromosome[rows], vcf.position[rows], vcf.ref_allele[rows],
	vcf.alt_allele[rows], vcf.gene[rows], vcf.effect[rows], vcf.notes[rows],
	vcf.alt[rows, cols], vcf.total[rows, cols], vcf.star[rows, cols], vcf.mapq[rows, cols], vcf.sidedness[rows, cols])
Base.getindex(vcf::VCF, rows, col::Int) = vcf[rows, [col]]
Base.getindex(vcf::VCF, row::Int, cols) = vcf[row, cols]
Base.size(vcf::VCF, args...) = Base.size(vcf.alt, args...)

mutable struct MinimalVCF
	sample::Vector{String}
	chromosome::Vector{String}
	position::Vector{Int32}
	ref_allele::Vector{String}
	alt_allele::Vector{String}
	alt::Matrix{Int32}
	total::Matrix{Int32}
	star::BitMatrix
end
MinimalVCF() = MinimalVCF([], [], [], [], [], zeros(0, 0), zeros(0, 0), falses(0, 0))
Base.getindex(vcf::MinimalVCF, rows, cols) = MinimalVCF(vcf.sample[cols],
	vcf.chromosome[rows], vcf.position[rows], vcf.ref_allele[rows],
	vcf.alt_allele[rows], vcf.alt[rows, cols], vcf.total[rows, cols],
	vcf.star[rows, cols])
Base.size(vcf::MinimalVCF, args...) = Base.size(vcf.alt, args...)

is_protein_altering(effect) = 
	r"Missense|[Ff]rameshift|Stopgain|Stoploss|Splice|Startloss|Promoter" in effect
is_missense(effect) = "issense" in effect
is_truncating(effect) = r"Frameshift|Stopgain|Stoploss|Splice|Startloss" in effect
is_coding_region(effect) = r"Synonymous|Missense|[Ff]rameshift|Stopgain|Stoploss|Splice|Startloss" in effect
is_indel(ref::AbstractString, alt::AbstractString) = length(ref) != length(alt)

function read_vcf_sample_names(vcf::Union{Cmd,AbstractString})
	vcf = open(vcf)
	for line in eachline(vcf)
		if startswith(line, '#'); continue; end
		headers = split(line, '\t')
		return headers[findone(headers, "NOTES")+1:end]
	end
	error("Header line was not found in VCF.")
end

function read_minimal_vcf(vcf_path::AbstractString)
	V = open(countlines, vcf_path)

	vcf = open(vcf_path)
	headers = split(readline(vcf), '\t')
	notes_col = findone(headers, "NOTES")
	samples = headers[notes_col+1:end]
	S = length(samples)

	chromosome = fill("", V)
	position = zeros(Int32, V)
	ref_allele = fill("", V)
	alt_allele = fill("", V)
	alt = zeros(Int32, V, S)
	total = zeros(Int32, V, S)
	star = falses(V, S)

	v = 0
	for line in eachline(vcf)
		cols = split(line, '\t')
		v += 1
		chromosome[v] = cols[1]
		position[v] = parse(Int32, cols[2])
		ref_allele[v] = cols[3]
		alt_allele[v] = cols[4]
		for s in 1:S
			parts = split(cols[notes_col + s], ':')		
			alt[v, s] = parse(Int32, parts[1])
			total[v, s] = parse(Int32, parts[2])
			star[v, s] = parts[end] == "*"
		end
	end
	return MinimalVCF(samples, chromosome, position, ref_allele, alt_allele, alt, total, star)
end

function read_vcf(vcf::Union{IO,AbstractString}; notes=false)
	d = readdlm(vcf, '\t')

	# Process header
	headers = d[1, :][:]; d = d[2:end, :];
	notes_col = findone(headers, "NOTES")
	samples = map(x -> "$x", headers[notes_col+1:end])
	filter!(s -> s != "BACKGROUND%", samples)
	S = length(samples)
	chromosome = String.(d[:, 1])
	position = convert(Vector{Int32}, d[:, 2])
	ref_allele = String.(d[:, 3])
	alt_allele = String.(d[:, 4])
	V = length(chromosome)

	gene_col = findone(headers, "GENE")
	gene = gene_col != nothing ? d[:, gene_col] : fill("", V)

	effect_col = findone(headers, "EFFECT")
	effect = effect_col != nothing ? d[:, effect_col] : fill("", V)

	alt = zeros(Int32, size(d, 1), S)
	total = zeros(Int32, size(d, 1), S)
	mapq = zeros(Int32, size(d, 1), S)
	sidedness = zeros(Int32, size(d, 1), S)
	star = falses(size(d, 1), S)

	for r in 1:size(d, 1)
		for s in 1:S
			parts = split(d[r, notes_col + s], ':')			
			if length(parts) <= 3
				# Call1
				star[r, s] = parts[end] == "*"
				alt[r, s] = parse(Int32, parts[1])
				total[r, s] = parse(Int32, parts[2])
			else 
				# Call2
				star[r, s] = parts[end] == "*"
				alt[r, s] = parse(Int32, parts[1])
				total[r, s] = parse(Int32, parts[2])
				mapq[r, s] = parse(Int32, parts[3])
				sidedness[r, s] = parse(Int32, parts[5])
			end				
		end
	end
	return VCF(samples, chromosome, position, ref_allele, alt_allele, gene, effect, notes ? d[:, notes_col] : fill("", V), alt, total, star, mapq, sidedness)
end

read_vcf(vcf_cmd::Base.AbstractCmd; kwargs...) = read_vcf(open(vcf_cmd); kwargs...)

# This function calculates a ctDNA fraction estimate based on the allele
# fraction of a single truncal mutation
function truncal_maf_to_cancer_frac_assuming_loh(alt, total, germline_ploidy)

	# The observed mutant allele fraction is drawn from a binomial
	# distribution Bin(n, p) where p = true mutant allele fraction.
	# Here we conservatively assume that the observed MAF may be a 95%
	# quantile anomaly.
	possible_mafs = 0:0.001:1
	maf = possible_mafs[searchsortedfirst(possible_mafs, true, by=(p -> cdf(Binomial(total, p), alt) < 0.95))]

	# Since we cannot detect copy number changes when cancer-%
	# is low, we conservatively always assume LOH for autosomes.
	cancer_frac = (germline_ploidy == 1) ? maf : 2 / (1 / maf + 1)
	return cancer_frac
end



function gene_mutated(vcf, gene; only_truncating=false, cohort=nothing, debug=false)
	V = size(vcf, 1)
	if cohort == nothing; cohort = trues(V); end
	if gene == "TERT promoter"
		rows = [vcf.chromosome[r] == "chr5" && 1_295_113 <= vcf.position[r] <= 1_295_135 for r in 1:V]
	else
		rows = [vcf.gene[r] == gene && is_protein_altering(vcf.effect[r]) && !("non-canonical isoform" in vcf.effect[r]) for r in 1:V]
	end

	mutated = falses(length(vcf.sample))
	for r in findall(rows); mutated .|= vcf.star[r, :]; end

	if debug
		for s in findall(mutated & cohort), r in findall(rows & vcf.star[:, s])
			@printf("%s: %s %s (%s:%d:%s>%s), %d / %d (%.1f%%)\n", vcf.sample[s], vcf.gene[r], vcf.effect[r], vcf.chromosome[r], vcf.position[r], vcf.ref_allele[r], vcf.alt_allele[r], vcf.alt[r, s], vcf.total[r, s], vcf.alt[r, s] / vcf.total[r, s] * 100)
		end
		@printf("--------------------\n")
	end

	return mutated[cohort]
end




function mutation_table(vcf; cohort=nothing, only_protein_altering=false, anonymous=nothing, output="mutations.tsv")
	if isfile(output); error("Output file '$output' already exists."); end
	if cohort == nothing; cohort = trues(length(vcf.sample)); end

	# Anonymous sample IDs can be provided as a vector or dictionary
	if anonymous == nothing
		anonymous = fill("", length(vcf.sample))
	elseif anonymous isa Dict
		anonymous = [get(anonymous, s, "") for s in vcf.sample]
	elseif anonymous isa AbstractVector
		@assert length(anonymous) == length(vcf.sample)
	end

	out = open(output, "w")
	@printf(out, "Sample\tAnon ID\tChrom\tPosition\tRef\tAlt\tGene\tEffect\tMutant allele fraction (mutant / total reads)\tBackground error rate\n")
	for r in 1:size(vcf, 1)
		if only_protein_altering && !is_protein_altering(vcf.effect[r])
			continue
		end
		
		for s in findall(cohort & vcf.star[r, :])
			@printf(out, "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s", vcf.sample[s],
				anonymous[s], vcf.chromosome[r], vcf.position[r],
				vcf.ref_allele[r], vcf.alt_allele[r], vcf.gene[r],
				vcf.effect[r])
			@printf(out, "\t%.1f%% (%d / %d)", vcf.alt[r, s] / vcf.total[r, s] * 100,
				vcf.alt[r, s], vcf.total[r, s])
			@printf(out, "\n")
		end
	end
	close(out)
end


end
