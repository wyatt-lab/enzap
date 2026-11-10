
module Plot2

export PyPlot, Matplotlib
export figure, RGB, hex, lighter_color
export line_plot, scatter_plot, bar_plot, violin_plot, area_plot, matrix_plot
export beeswarm_plot, stacked_bar_plot, box_plot, heatmap
export rectangle, arc, ellipse, errorbar, text
export matrix_plot
export xlim, ylim, xticks, yticks, subplot, ylabel, xlabel, title, grid
export genome_scatter_plot, genome_segments_plot, genome_labels, genome_plot_config
export suptitle, arrow_plot, histogram_plot
export get_ylim, get_xlim, square_lims, diagonal
export start_drawing, stop_drawing
export density_beeswarm
export grouped_barplot
using Helpers, Printf, DelimitedFiles, CopyNum, KernelDensity, Random
import PyCall

global Matplotlib = nothing;
global PyPlot = nothing;
global drawing_figure = false;
global hg38_chr_sizes = [
	"chr1"   248956422; "chr2"   242193529;	"chr3"   198295559;
	"chr4"   190214555;	"chr5"   181538259;	"chr6"   170805979;
	"chr7"   159345973;	"chr8"   145138636;	"chr9"   138394717;
	"chr10"  133797422;	"chr11"  135086622;	"chr12"  133275309;
	"chr13"  114364328;	"chr14"  107043718;	"chr15"  101991189;
	"chr16"   90338345;	"chr17"   83257441;	"chr18"   80373285;
	"chr19"   58617616;	"chr20"   64444167;	"chr21"   46709983;
	"chr22"   50818468;	"chrX"   156040895;	"chrY"    57227415;	
	"chrM"       16569
]

function figure(block::Function, path::String; size=(5, 3), dpi=300)
	global drawing_figure, Matplotlib, PyPlot;
	if drawing_figure == true
		block()
	else
		path = expanduser(path)
		if !isdir(dirname(abspath(path)))
			error("Cannot create figure $path. Directory does not exist.")
		end

		# Import matplotlib.pyplot if we have not done so yet
		if PyPlot == nothing
			Matplotlib = PyCall.pyimport("matplotlib")
			PyPlot = PyCall.pyimport("matplotlib.pyplot")
		end

		drawing_figure = true

		pyplot_fig = PyPlot.figure(figsize=size,dpi=dpi)
		try
			transparent_bg = r"\.(pdf|PDF|svg|SVG)$" in path
			block()
			PyPlot.tight_layout()
			PyPlot.savefig(path, transparent=transparent_bg)
			println(stderr, "Figure saved at $(path)")
		finally
			PyPlot.close()
			drawing_figure = false
		end
	end
	return nothing
end

figure(block::Function; kwargs...) = figure(block, expanduser("~/plot.pdf"); kwargs...)

xlabel(label::AbstractString) = PyPlot.xlabel(label)
ylabel(label::AbstractString) = PyPlot.ylabel(label)
function title(text::AbstractString; offset=0, location="center", size="large")
	PyPlot.title(text, pad=offset, loc=location, size=size)
end

function suptitle(text::AbstractString; fontsize=15)
	PyPlot.suptitle(text, fontsize=fontsize)
end

function xticks(ticks::AbstractVector; rotation=90, spine=true, 
				grid=nothing, labels=nothing, label_rotation=0)
	if all(x -> x isa AbstractString, ticks)
		PyPlot.xticks(1:length(ticks), ticks, rotation=rotation)
		PyPlot.tick_params(axis="x", length=0)
	else
		PyPlot.xticks(ticks)
	end
	PyPlot.gca().spines["bottom"].set_visible(spine)
	PyPlot.gca().spines["top"].set_visible(false)
	if grid == true
		PyPlot.grid(true, axis="x", linestyle="dashed")
	elseif grid == false
		PyPlot.grid(false, axis="x")
	end
	if labels != nothing
		PyPlot.gca().set_xticklabels(labels, rotation=label_rotation)
	end
end

function yticks(ticks::AbstractVector; spine=true, grid=nothing, labels=nothing)
	if all(x -> x isa AbstractString, ticks)
		PyPlot.yticks(1:length(ticks), ticks)
		PyPlot.tick_params(axis="y", length=0)
	else
		PyPlot.yticks(ticks)
	end
	PyPlot.gca().spines["left"].set_visible(spine)
	PyPlot.gca().spines["right"].set_visible(false)
	if grid == true
		PyPlot.grid(true, axis="y", linestyle="dashed")
	elseif grid == false
		PyPlot.grid(false, axis="y")
	end
	if labels != nothing
		PyPlot.gca().set_yticklabels(labels)
	end
end

#remove_xtick_labels() = PyPlot.gca().set_xticklabels([]); PyPlot.gca().set_xticklabels([], minor=true)
#remove_ytick_labels() = PyPlot.gca().set_yticklabels([]); PyPlot.gca().set_yticklabels([], minor=true)

grid(args...; kwargs...) = PyPlot.grid(args...; kwargs...)

function xlim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.xlim(left=low); end
	if isfinite(high); PyPlot.xlim(right=high); end
	if log > 1; PyPlot.semilogx(base=log); end
end

function ylim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.ylim(bottom=low); end
	if isfinite(high); PyPlot.ylim(top=high); end
	if log > 1; PyPlot.semilogy(base=log); end
end
get_ylim() = PyPlot.gca().get_ylim()
get_ylim(axis) = PyPlot.gca().get_ylim()[axis]
get_xlim() = PyPlot.gca().get_xlim()
get_xlim(axis) = PyPlot.gca().get_xlim()[axis]

function square_lims() # extend lims so that axes show the same range
    mn = min(get_xlim(1), get_ylim(1))
    mx = max(get_xlim(2), get_ylim(2))
    xlim(mn, mx)
    ylim(mn, mx)
end
function diagonal() # draw a diagonal grid line
    xmn, xmx = get_xlim()
    ymn, ymx = get_ylim()
    mn = max(xmn, ymn)
    mx = min(xmx, ymx)
    line_plot([mn, mx], [mn, mx], color=RGB(175), style="dashed", line_width=0.8)
    xlim(xmn, xmx)
    ylim(ymn, ymx)
end

function subplot(panel::Integer, rows::Integer)
	PyPlot.subplot(rows, 1, panel)
end

function subplot(panel::Integer; rows::Integer=1, cols::Integer=1)
	PyPlot.subplot(rows, cols, panel)
end

function subplot(panel::Tuple{Integer, Integer}; rows::Integer=1, cols::Integer=1, 
                 rowspan::Integer=1, colspan::Integer=1)
	PyPlot.subplot2grid((rows, cols), panel .- 1, rowspan=rowspan, colspan=colspan)
end

struct RGB; r::UInt8; g::UInt8; b::UInt8; end
RGB(gray::Integer) = RGB(gray, gray, gray)
Base.broadcastable(r::RGB) = Ref(r)
hex(c::RGB) = @sprintf("#%02X%02X%02X", c.r, c.g, c.b)
hex(colors::Vector{RGB}) = hex.(colors)
function lighter_color(rgb::RGB; factor=0.2)
	RGB(min( round(rgb.r+(255-rgb.r)*factor, digits=0), 255 ), 
		min( round(rgb.g+(255-rgb.g)*factor, digits=0), 255 ), 
		min( round(rgb.b+(255-rgb.b)*factor, digits=0), 255 ))
end

function text(x, y, str; size=10, h_align="center", v_align="center", color=RGB(0), rotation="horizontal")
	PyPlot.text(x, y, str, fontsize=size, horizontalalignment=h_align, verticalalignment=v_align, color=hex(color), rotation=rotation)
end

function line_plot(x::AbstractVector, y::AbstractVector; color=RGB(0,0,0), line_width=1, capstyle="butt", style="solid", zorder=10)
	figure() do
		PyPlot.plot(x, y, color=hex(color), linewidth=line_width, solid_capstyle=capstyle, linestyle=style, zorder=zorder)
		PyPlot.grid(true, axis="both", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function errorbar(args...; kwargs...)
	figure() do
		PyPlot.errorbar(args...; kwargs...)
	end
end

function genome_plot_config(chr_names, chr_lens)
	PyPlot.grid(true, axis="x", which="major", linestyle="dashed")
	PyPlot.grid(false, axis="y", which="both")   # No y-axis grid lines
	PyPlot.xlim(0, sum(chr_lens))
	PyPlot.tick_params(axis="x", which="both", length=0)
	PyPlot.gca().set_xticks(cumsum(chr_lens))
	PyPlot.gca().set_xticklabels([])
	PyPlot.gca().set_xticks(cumsum(chr_lens) .- chr_lens ./ 2, minor=true)
	PyPlot.gca().set_xticklabels(
		[r"(19|21|M)$" in chr ? "" : replace(chr, "chr", "")
		for chr in chr_names], minor=true)
end

function genome_scatter_plot(chromosome::AbstractVector,
	position::AbstractVector, value::AbstractVector;
	chr_sizes=hg38_chr_sizes, color=RGB(0), line=false,
	size=2, marker="o")

	@assert length(chromosome) == length(position) == length(value)
	if isa(color, AbstractArray)
		@assert length(color) == length(chromosome)
	end

	d = isa(chr_sizes, String) ? readdlm(expanduser(chr_sizes)) : chr_sizes
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])

	keep = [c in chr_names for c in chromosome]
	if sum(keep) != length(chromosome)
		println(stderr, "WARNING: All chromosomes not found in reference: ")
		println(join(unique(chromosome[!keep]), ", "))
	end
	chromosome = chromosome[keep]
	position = position[keep]
	value = value[keep]
	if isa(color, AbstractArray)
		color = color[keep]
	end

	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	gpos = [chr_start[chromosome[k]] + position[k]
			for k in 1:length(chromosome)]

	figure() do
		if line
			order = sortperm(gpos)
			area_plot(gpos[order], value[order], color=color)
		else
			scatter_plot(gpos, value, color=color, size=size, marker=marker)
		end

		genome_plot_config(chr_names, chr_lens)
	end
end

function genome_segments_plot(segments::Vector{Segment},
	seg_values::Vector; chr_sizes=hg38_chr_sizes,
	color=RGB(255, 0, 0), line_width=2)

	@assert(length(segments) == length(seg_values))
	d = isa(chr_sizes, String) ? readdlm(expanduser(chr_sizes)) : chr_sizes
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	figure() do
		for (i, (seg, value)) in enumerate(zip(segments, seg_values))
			start = chr_start[seg.chromosome] + seg.start
			stop = chr_start[seg.chromosome] + seg.stop
			c = typeof(color) <: AbstractVector ? color[i] : color
			line_plot([start, stop], fill(value, 2),
					color=c, line_width=line_width)
		end
		genome_plot_config(chr_names, chr_lens)
	end
end

function genome_labels(chromosome::AbstractVector,
	position::AbstractVector, ypos::AbstractVector, label::AbstractVector;
	chr_sizes=hg38_chr_sizes, color=RGB(0), size=8)

	@assert(length(chromosome) == length(position) == length(ypos) == length(label))

	d = isa(chr_sizes, String) ? readdlm(expanduser(chr_sizes)) : chr_sizes
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
		for c in 1:length(chr_names))

	gpos = [chr_start[chromosome[k]] + position[k]
		for k in 1:length(chromosome)]

	figure() do
		for k in 1:length(label)
			text(gpos[k], ypos[k], label[k]; size=size, v_align="top", rotation="vertical", color=color)
			#text(x, y, str; size=10, h_align="center", v_align="center", color=RGB(0), rotation="horizontal")
		end
	end

end

# using Distributions
# function genome_kde_plot(chromosome::AbstractVector, position::AbstractVector;
# 	chr_sizes=hg38_chr_sizes, sigma=1e6, color=RGB(0))

# 	@assert(length(chromosome) == length(position))
# 	d = isa(chr_sizes, String) ? readdlm(expanduser(chr_sizes)) : chr_sizes
# 	chr_names = d[:, 1]
# 	chr_lens = Int32.(d[:, 2])
# 	chr_start = map(c -> sum(chr_lens[1:c-1]), eachindex(chr_lens))

# 	gpos = map(1:length(chromosome)) do k
# 		chr = findone(chr_names, chromosome[k])
# 		start = chr_start[chr]
# 		pos = start + position[k]
# 		stop = start + chr_lens[chr]
# 		return [pos, start, stop]
# 	end
	
# 	x = 0:1e6:(chr_start[end]+chr_lens[end])
# 	y = zeros(length(x))
# 	for gp in gpos
# 		i = gp[2] .< x .< gp[3]
# 		yy = pdf(Normal(gp[1], sigma), x[i])
# 		y[i] += yy
# 	end

# 	figure() do
# 		line_plot(x, y, color=color)
# 		genome_plot_config(chr_names, chr_lens)
# 	end
# end

# function genome_plot_annotation(chromosome::AbstractVector, position::AbstractVector, annotation::AbstractVector; 
# 	chr_sizes=hg38_chr_sizes, color=RGB(255,0,0), fontsize=10, line_width=1, rotation="horizontal")

# 	@assert(length(chromosome) == length(position))
# 	d = isa(chr_sizes, String) ? readdlm(expanduser(chr_sizes)) : chr_sizes
# 	chr_names = d[:, 1]
# 	chr_lens = Int32.(d[:, 2])
# 	chr_start = Dict(chr_names[c] => sum(chr_lens[1:c-1])
# 		for c in 1:length(chr_names))

# 	gpos = [chr_start[chromosome[k]] + position[k]
# 			for k in 1:length(chromosome)]
	
# 	ylim = [l for l in PyPlot.gca().get_ylim()]

# 	for (i, pos) in enumerate(gpos)
# 		line_plot([pos, pos], ylim, color=color, line_width=line_width)
# 		text(pos, ylim[2], annotation[i], size=fontsize, v_align="bottom", rotation=rotation)
# 	end
# end

using KernelDensity
function kde_plot(x::AbstractArray; bandwidth=0.5, line_width=1, npoints=2048, color=RGB(0))
	U = kde(Float64.(x), bandwidth=bandwidth, npoints=npoints)
	line_plot(U.x, U.density, line_width=line_width, color=color)
end

function scatter_plot(x::AbstractVector, y::AbstractVector;
	size=2, linewidths=1.0, color=RGB(0,0,0), opacity=1.0, border_color=nothing, fill_color=nothing, marker="o", 
	xerror=nothing, yerror=nothing, error_line_width=1, error_color=RGB(0,0,0), grid=true, zorder=1,rasterized=false)

	@assert(length(x) == length(y))

	# This is a workaround for the fact that Matplotlib misindexes the
	# "error_color" vector if some error ranges are NaN.
	if xerror != nothing
		@assert(xerror isa AbstractMatrix && Base.size(xerror, 1) == length(x))
		for k in 1:Base.size(xerror, 1)
			if any(isnan, xerror[k, :]); xerror[k, :] .= 0; end
		end
	end
	if yerror != nothing
		@assert(yerror isa AbstractMatrix && Base.size(yerror, 1) == length(x))
		for k in 1:Base.size(yerror, 1)
			if any(isnan, yerror[k, :]); yerror[k, :] .= 0; end
		end
	end

	figure() do
		if xerror isa AbstractMatrix
			PyPlot.hlines(y, xmin=xerror[:, 1], xmax=xerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end
		if yerror isa AbstractMatrix
			PyPlot.vlines(x, ymin=yerror[:, 1], ymax=yerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end

		if grid == true
			PyPlot.grid(true, axis="both", linestyle="dashed")
				#color=hex.(RGB(220)), linewidth=0.1)
		else
			PyPlot.grid(false)
		end
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		edgecolors = all(isa.(border_color, RGB)) ? hex(border_color) : "none"
		facecolors = isnothing(fill_color) ? hex(color) : isa(fill_color, RGB) ? hex(fill_color) : fill_color
		color = facecolors=="none" ? "none" : hex(color)
		PyPlot.scatter(x, y, s=size, linewidths=linewidths, marker=marker, c=color,
			alpha=opacity, edgecolors=edgecolors, facecolors=facecolors, zorder=zorder,rasterized=rasterized)
	end
end

function bar_plot(values::AbstractVector; box_width=0.8, color=RGB(0,0,0),opacity=1.0)
	figure() do
		PyPlot.bar(1:length(values), values, color=hex(color),
			width=box_width, linewidth=0,
			bottom=0.0000000000000001, zorder=10,alpha=opacity)
		PyPlot.xlim(0, length(values) + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function stacked_bar_plot(values::AbstractMatrix; box_width=0.8, colors=[])
	B = size(values, 1)    # Number of bars
	L = size(values, 2)    # Number of levels

	if isempty(colors)
		colors = [RGB(round(Int, (l - 1) / L) * 255) for l in 1:L]
	end

	@assert(colors isa AbstractVector{RGB})
	if length(colors) != L
		error("Bar plot has $L levels but $(length(colors)) colors were provided in 'colors' keyword argument.")
	end

	figure() do
		bottoms = zeros(B)
		PyPlot.xlim(0, B + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for level in 1:L
			PyPlot.bar(1:B, values[:, level], width=box_width, bottom=bottoms, color=hex(colors[level]), zorder=10)
			bottoms .+= values[:, level]
		end
	end
end

function box_plot(groups...; dots="all", box_width=0.8, dot_size=20, beeswarm_edge=true, color=RGB(200), vertical=true)
	figure() do
		if !vertical; @assert dots != "all"; end
		boxes = PyPlot.boxplot(groups, showfliers=false, widths=box_width,
			patch_artist=true, vert=vertical)
		if isa(color, AbstractVector)
			@assert length(groups) == length(color)
			for b in 1:length(groups)
				boxes["boxes"][b].set_facecolor(hex(color[b]))
			end
		else
			for box in boxes["boxes"]
				box.set_facecolor(hex(color))
			end
		end
		for med in boxes["medians"]; med.set_color("black"); end
		if dots == "all"
			beeswarm_plot(groups...; point_size=dot_size, marker="o", color=RGB(255), edgecolor=beeswarm_edge ? "black" : "none")
		end
	end
end

function area_plot(x::AbstractVector, y::AbstractVector; xbottom=nothing, ybottom=nothing, 
			color=RGB(0), alpha=1.0, zorder=1)

	figure() do
		if xbottom == nothing && ybottom == nothing
			PyPlot.fill_between(x, y, facecolor=hex(color), alpha=alpha, zorder=zorder)
		end
		if xbottom isa Number || xbottom isa AbstractVector
			PyPlot.fill_betweenx(y, xbottom, x, facecolor=hex(color), alpha=alpha, zorder=zorder)
		end
		if ybottom isa Number || ybottom isa AbstractVector
			PyPlot.fill_between(x, ybottom, y, facecolor=hex(color), alpha=alpha, zorder=zorder)
		end
	end
end

function histogram_plot(x; bins=nothing, nbins=nothing, bwidth=nothing, pmf=false, kwargs...)

	if bins !== nothing
		@assert nbins === nothing && bwidth === nothing
	else
		@assert !(nbins !== nothing && bwidth !== nothing)
	end

	x = filter(isfinite, x)
	if isempty(x)
		if bins !== nothing
			return bins, zeros(eltype(x), length(bins)-1)
		elseif nbins !== nothing
			return nothing, zeros(eltype(x), nbins)
		else
			return nothing, nothing
		end
	end

	if bins === nothing
		mn, mx = extrema(x)
		if bwidth !== nothing
			bins = range(mn, mx, step=bwidth)
		else
			nbins === nothing && (nbins = 100)
			bins = range(mn, mx, length=nbins+1)
		end
	end

	bincounts = zeros(length(bins)-1)
	for xx in x
		if xx < bins[1]; continue; end
		i = searchsortedfirst(bins, xx)-1
		i==0 && (i = 1)
		if i > length(bincounts); continue; end
		bincounts[i] += 1
	end
					
	if pmf; bincounts ./= sum(bincounts); end
	bins = [(bins[i] + bins[i+1]) / 2 for i in 1:length(bins)-1]
	line_plot(bins, bincounts; kwargs...)
end

function violin_plot(groups...; color=[], halved=false, kwargs...)
	lowest = minimum(minimum.(groups))
	highest = maximum(maximum.(groups))
	lowest -= 0.5 * (highest - lowest) 
	highest += 0.5 * (highest - lowest)
	grid = range(lowest, highest, length=500)
	density = zeros(length(grid), length(groups))
	for g in 1:length(groups)
		density[:, g] = kde(groups[g], grid; kwargs...).density
	end
	violin_plot(grid, density; color=color, halved=halved)
end

function violin_plot(bins::AbstractVector, density::Matrix; color=[], halved=false)

	V = size(density, 2);   # How many violins?

	bins = Float64.(bins); density = Float64.(density)
	for v in 1:V
		density[:, v] /= maximum(density[:, v]) * 2.05
	end

	if color isa RGB
		color = fill(color, V)
	elseif isempty(color)
		color = map(v -> RGB(200), 1:V)
	elseif color isa AbstractVector
		@assert(length(color) == V)
	end

	figure() do
		PyPlot.xlim(0.25, V + 0.75)
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for v in 1:V
			if halved
				x = vcat(v .- density[:, v], v, v, v .- density[1, v])
				y = vcat(bins, bins[end], bins[1], bins[1])
			else
				x = vcat(v .- density[:, v], v .+ density[end:-1:1, v])
				y = vcat(bins, bins[end:-1:1])
			end
			PyPlot.fill(x, y, hex(color[v]), zorder=10)
		end
	end
end

function rectangle(x1::Real, x2::Real, y1::Real, y2::Real; 
			edge_width=0, color=RGB(0, 0, 0), opacity=1.0, zorder=1)
	width = abs(x2 - x1); height = abs(y2 - y1);
	x = min(x1, x2); y = min(y1, y2);
	rect = PyPlot.matplotlib.patches.Rectangle((x, y), width, height, linewidth=edge_width, 
					facecolor=hex(color), alpha=opacity, zorder=zorder)
	PyPlot.gca().add_patch(rect)
end

function arc(x1::Real, x2::Real, y1::Real, y2::Real; color=RGB(0))
	@assert(y2 >= y1)
	patch = PyPlot.matplotlib.patches.Arc(((x1 + x2) / 2, y1), abs(x2 - x1), y2 - y1, theta2=180)
	PyPlot.gca().add_patch(patch)
end

function ellipse(x::Real, y::Real, rx::Real, ry::Real; color=RGB(0, 0, 0))
	patch = PyPlot.matplotlib.patches.Ellipse([x, y], 2*rx, 2*ry,
		facecolor=hex(color))
	PyPlot.gca().add_patch(patch)
end

function beeswarm_plot(groups...; point_size=10, color=RGB(0,0,0), marker=".", edgecolor="none")
	x = zeros(0); y = zeros(0); face_colors = RGB[];
	for (group, values) in enumerate(groups)
		append!(x, group .+ clamp.(randn(length(values)) / 10, -0.4, 0.4))
		append!(y, values)
		if color isa AbstractVector
			append!(face_colors, color[group])
		else
			append!(face_colors, fill(color, length(values)))
		end
	end
	figure() do
		xticks([])
		#PyPlot.tick_params(axis="x", length=0)
		#PyPlot.grid(false)
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
		PyPlot.scatter(x, y, s=point_size, c=hex.(face_colors), marker=marker,
			edgecolors=edgecolor, zorder=10)
	end
end

function density_jitter(x, GAP;seed=1234)
	Random.seed!(seed)
    #shuffle to avoid "worms"
    rand_ord = randperm(length(x))
    x = x .- minimum(x) .+ GAP
    spans = [(round(Int, 100_000 * (a - GAP)) + 1):(round(Int, 100_000 * (a + GAP)) + 1)  for a in x[rand_ord]]
    span_size = round(Int, 100_000 * (maximum(x) + GAP)) + 1

	levels = fill(NaN, length(x))
	for level in 0:1000
		occupied = falses(span_size)
        neg_occupied = falses(span_size)
		for r in 1:length(x)
			if isfinite(levels[r]); continue; end       # Already has a level
			if any(occupied[spans[r]]) # Does not fit
    
                if level .== 0;continue;end
                if any(neg_occupied[spans[r]])# Does not fit
                    continue
                else
                    levels[r] = -1*level	# This fits
                    neg_occupied[spans[r]] .= true
                end
            else
                levels[r] = level	# This fits
                occupied[spans[r]] .= true
            end
		end

		if all(isfinite, levels); break; end
	end
    out = fill(0.0, length(x))
    for (i, k) in enumerate(rand_ord)
        out[k] = levels[i]
    end
    out ./= max(1,maximum(abs.(out)))   # Scale between -1 and 1
    return out
end



function density_beeswarm(groups...;dot_size="auto",color="black",edgecolor="black",marker="o",seed=1234,median_line=true,data_range=NaN)
    fig = PyPlot.gcf()

    if dot_size == "auto"
        estimates=[]
        for (i,values) in enumerate(groups)
			if isempty(values)
				push!(estimates,Inf)
				continue
			end
            kd=kde(Float64.(values))
            peak=findfirst(kd.density.==maximum(kd.density))
            window=abs(maximum(values)-minimum(values))*0.02
            n=max(1,sum(kd.x[peak]-window .< values .< kd.x[peak] + window))
            push!(estimates,1/sqrt(n)*0.05)
        end
        dot_size=minimum(estimates)
        print("Estimated dot size $(dot_size)")
	end
      
	if isnan(data_range)
		max_val=maximum(vcat(groups...))
		min_val=minimum(vcat(groups...))
    	gap=(max_val-min_val)* dot_size / 2
	else
		gap=data_range* dot_size / 2
	end
    y=[]
    x=[]
    max_levels=[]
    for (i,values) in enumerate(groups)
		if isempty(values)
			push!(x,[])
			push!(y,[])
			push!(max_levels,0)
			continue
		end
        jitter=density_jitter(values,gap,seed=seed)
        max_level=ceil((length(unique(jitter)) - 1)/2)
        push!(x,jitter)
        push!(y,values)
        push!(max_levels,max_level)
    end
    scales=max_levels ./ maximum(max_levels)
    y=vcat(y...)
    
    ax = PyPlot.subplot(1,1,1) 
    #set the axis limits to scale the dot size
    PyPlot.ylim(minimum(y)-0.02*abs(maximum(y)-minimum(y)),maximum(y)+0.02*abs(maximum(y)-minimum(y)))
    PyPlot.xlim(0.5,length(groups)+0.5)
    #point size is scaled between maximum and minimum y values
    s=(0.5 *dot_size * (ax.transData.transform((0,maximum(y)))[2]-ax.transData.transform((0,minimum(y)))[2])) ^2
    
    #scale the X-axis gap
    x_gab=2*gap*(ax.get_xlim()[2]-ax.get_xlim()[1]) / (ax.get_ylim()[2]-ax.get_ylim()[1])
    #scale X-axis gap by figure width
    x_gab *= (fig.get_size_inches()[2] / fig.get_size_inches()[1])
    width=min(0.45,x_gab*maximum(max_levels))
    x=vcat([x[i] .* width .* scales[i] .+ i for i in 1:length(max_levels)]...)
    scatter_plot = ax.scatter(x, y,s=s,c=color,marker=marker,edgecolors=edgecolor, zorder=15)
    ax.set_xticks(1:length(groups))

	if median_line
        y_pos=[median(values) for values in groups]
        xmin=[i-0.25 for i in 1:length(groups)]
        xmax=[i+0.25 for i in 1:length(groups)]
		PyPlot.hlines(y_pos,xmin=xmin,xmax=xmax,linewidth=3,zorder=16)

	end

end

#added by Tuomo
function grouped_barplot(data;labels=[],colors=[],bar_width="auto")
    if isempty(labels)
        labels=1:length(data[1])
    end
    figure() do
        #bar_width = 0.4
        x_pos = 1:length(labels)
		if bar_width == "auto";bar_width=1/length(data)*0.8;end
        # Loop through each group of bars and create a set of bars for each group
        for (i,dt) in enumerate(data)
            # Calculate the x-axis position for the current group of bars
            pos = x_pos .+ (i * bar_width)
            # Create a set of bars for the current group
            if !isempty(colors)
                PyPlot.bar(pos, dt, width=bar_width, label=labels,color=hex(colors[i]))
            else
                PyPlot.bar(pos, dt, width=bar_width, label=labels)
            end
        end
        PyPlot.xticks(x_pos .+ ((length(data)) / 2) * bar_width .+ bar_width/2 ,labels=labels,rotation=45,ha="right")
    end
end


function arrow_plot(x, y, dx, dy; width=0.001, 
	head_width=3*width, head_length=1.5*head_width,
	color=RGB(0), zorder=2, length_includes_head=false)
	PyPlot.arrow(x, y, dx, dy, width=4e6, head_width=head_width, 
		head_length=head_length, color=Plot2.hex(color), zorder=zorder, 
		length_includes_head=length_includes_head)
end

function matrix_plot(colors::Matrix{RGB}; path="~/plot.svg", cell_width=50, cell_height=50, merge=false)
	R, C = size(colors)
	canvas_w = C * cell_width
	canvas_h = R * cell_height
	svg = open(expanduser(path), "w")
	write(svg, """<svg width="$canvas_w" height="$canvas_h">""")
	painted = falses(size(colors))
	for r in 1:R, c in 1:C
		if painted[r, c]; continue; end

		color = colors[r, c]

		# Check if we can merge same-colored adjacent cells
		rowspan = 1; colspan = 1
		if merge
			while r + rowspan <= size(colors, 1) && colors[r+rowspan, c] == color; rowspan += 1; end
			while c + colspan <= size(colors, 2) && colors[r, c+colspan] == color; colspan += 1; end

			# Currently we only ever merge along rows or columns, but not both
			if rowspan > colspan
				colspan = 1
			else
				rowspan = 1
			end
		end

		x = (c - 1) * cell_width
		y = (r - 1) * cell_height
		rgb = "rgb($(color.r), $(color.g), $(color.b))"
		write(svg, """<rect x="$x" y="$y" width="$(cell_width * colspan)" height="$(cell_height * rowspan)" style="stroke: none; fill: $rgb" />""")

		painted[r:r+rowspan-1, c:c+colspan-1] .= true
	end
	write(svg, "</svg>")
	close(svg)
end

function heatmap(values::AbstractMatrix; colormap="viridis")
	figure() do
		PyPlot.imshow(values; cmap=colormap, aspect="auto")
	end
end

#function heatmap(colors::Matrix{RGB}; path="~/plot.png")
#	path = expanduser(path)
#	rgb = zeros(UInt8, size(colors, 1), size(colors, 2), 3)
#	rgb[:, :, 1] = map(c -> c.r, colors)
#	rgb[:, :, 2] = map(c -> c.g, colors)
#	rgb[:, :, 3] = map(c -> c.b, colors)
#	PNGFiles.save(path, rgb)
#end



# Functions needed for joint usage of Plot2 and the public package PyPlot
function start_drawing()
	global drawing_figure, PyPlot
	if PyPlot == nothing
		PyPlot = PyCall.pyimport("matplotlib.pyplot")
	end
	drawing_figure = true
end
stop_drawing() = global drawing_figure = false



############
# COLORS #
############

export GROUP_PALETTE, diverging_color, diverging_color2, color_scale, circular_color, hsl_to_rgb, hsl_color_palette;
const GROUP_PALETTE = [RGB(255, 0, 0), RGB(0, 0, 255), RGB(0, 204, 0),
	RGB(255, 0, 255), RGB(255, 153, 51), RGB(0, 230, 230), RGB(230, 230, 0),
	RGB(51, 153, 255), RGB(174, 230, 100), RGB(255, 0, 127), RGB(255, 153, 153),
	RGB(153, 51, 255), RGB(184, 134, 11), RGB(128, 128, 200),
	RGB(238, 130, 238), RGB(75, 0, 130)];

function hsl_to_rgb(h, s, l)
	@assert 0 <= s <= 1
	@assert 0 <= l <= 1
	function hue_to_rgb(p, q, t)
		if t < 0
			t += 1
		elseif t > 1
			t -= 1
		end
		if t < 1/6
			return p + (q - p) * 6 * t
		elseif t < 1/2
			return q
		elseif t < 2/3
			return p + (q - p) * (2/3 - t) * 6
		else
			return p
		end
	end

	q = l < 0.5 ? l * (1 + s) : l + s - l * s
	p = 2 * l - q
	h = h / 360
	r = hue_to_rgb(p, q, h + 1/3)
	g = hue_to_rgb(p, q, h)
	b = hue_to_rgb(p, q, h - 1/3)
	
	return RGB(round.(Int, [r * 255, g * 255, b * 255])...)
end
function hsl_color_palette(n, h0, s, l)
	@assert 0 <= h0 <= 360 
	h = collect(range(0, 360, length=n+1)[1:end-1])
	for i in 1:n
		h[i] += h0
		if h[i] >= 360
			h[i] %= 360
		end
	end
	return hsl_to_rgb.(h, s, l)
end

# -1=minus, 0=black, 1=plus
function diverging_color(minus::RGB, plus::RGB, val)
	if isnan(val)
		RGB(220)
	elseif val < 0
		RGB(round(Int, -val*minus.r), round(Int, -val*minus.g), round(Int, -val*minus.b))
	else
		RGB(round(Int, val*plus.r), round(Int, val*plus.g), round(Int, val*plus.b))
	end
end
function diverging_color(min::RGB, max::RGB, vec::AbstractVector; top_q=1.0)
    scaled = vec .- median(filter(isfinite, vec))
    scaled ./= quantile(abs.(filter(isfinite, scaled)), top_q)
	scaled = clamp.(scaled, -1, 1)
    return diverging_color.(min, max, scaled);
end
# linear interpolation of colors between values
function color_scale(x::T where T <: Number, values, colors; nacol=RGB(255))
    @assert length(values) == length(colors)
    if !isfinite(x); return nacol; end
    first = findlast(x .>= values)
    if isnothing(first); return colors[1]; end
    if first == length(values); return colors[end]; end
    col1, col2 = colors[first:first+1]
    val1, val2 = values[first:first+1]
    sc = (x-val1) / (val2-val1) - 1e-9
    return col2*sc + col1*(1-sc)
end
function color_scale(x::T where T <: AbstractArray, values, colors; nacol=RGB(255))
    return [color_scale(xx, values, colors, nacol=nacol) for xx=x]
end

# -1=minus, 1=plus
using Statistics
function diverging_color2(minus::RGB, plus::RGB, val; nan=RGB(220))
	if isnan(val); return nan; end
	f = (val + 1) / 2
	RGB(round(Int, f*plus.r + (1-f)*minus.r), 
		round(Int, f*plus.g + (1-f)*minus.g), 
		round(Int, f*plus.b + (1-f)*minus.b))
end
function diverging_color2(min::RGB, max::RGB, vec::AbstractVector; top_q=NaN, nan=RGB(220))
	if isnan(top_q)
		scaled = vec .- minimum(filter(isfinite, vec))
		scaled ./= maximum(filter(isfinite, scaled))
		scaled .*= 2
		scaled .-= 1
	else
		scaled = vec .- mean(filter(isfinite, vec))
		scaled ./= quantile(abs.(filter(isfinite, scaled)), top_q)
	end
    return diverging_color2.(min, max, clamp.(scaled, -1, 1); nan=nan);
end

# Create a circular color scale with <colors>, and find color of <t> ∈ [0, 1]
function circular_color(colors::Vector{RGB}, t; nacol=RGB(255))
	isnan(t) && return nacol
    @assert 0 <= t <= 1
    n = length(colors)
    t *= n
    prev = floor(Int, t) + 1
    next = (prev % n) + 1
    f = t % 1
    return (1-f) * colors[prev] + f * colors[next]
end

function density_color(a, b; bandwidth=(1.0,1.0), top_q=1.0,
        color=[RGB(0), RGB(0,0,255), RGB(0,165,255), RGB(0,255,165), 
        RGB(0,255,0), RGB(165,255,0), RGB(255,165,0), RGB(255,0,0)])
    
    valid = isfinite.(a) & isfinite.(b)

    a_bw = std(a[valid]) * bandwidth[1] / 3
    b_bw = std(b[valid]) * bandwidth[2] / 3

    ik = InterpKDE(kde((a[valid], b[valid]), bandwidth=(a_bw, b_bw)))
    dens = [pdf(ik, ee, cc) for (ee, cc)=zip(a, b)]

    mx = quantile(filter(isfinite, dens), top_q)
    mn = quantile(filter(isfinite, dens), 1-top_q)
    vals = mn .+ ((mx-mn)/length(color)) .* (0:length(color)-1)
    return color_scale(dens, vals, color)
end

import Base.*, Base.+
*(v::AbstractFloat, c::RGB) = RGB(round(Int, v*c.r), round(Int, v*c.g), round(Int, v*c.b))
*(c::RGB, v::AbstractFloat) = RGB(round(Int, v*c.r), round(Int, v*c.g), round(Int, v*c.b))
+(a::RGB, b::RGB) = RGB(Int(a.r)+Int(b.r), Int(a.g)+Int(b.g), Int(a.b)+Int(b.b)) # convert to int so if we would overflow, we give error instead
Statistics.mean(a::RGB, b::RGB) = 0.5*a + 0.5*b

end

