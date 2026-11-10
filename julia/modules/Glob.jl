
module Glob

export glob

function glob(pattern::String)
	if occursin(r"[\^\$\\]", pattern)
		error("Special symbols ^ \$ \\ are not supported in glob pattern.")
	end
	pattern = expanduser(pattern)
	dir = dirname(pattern)
	if "*" in dir; error("* is not supported for directories yet"); end
	files = readdir(dir == "" ? "." : dir)

	rx = Regex("^$(replace(basename(pattern), "*" => ".*", "." => "\\.", "?" => "."))\$")
	return [dir == "" ? f : "$dir/$f" for f in files if occursin(rx, f)]
end

glob(patterns...) = vcat(glob.(patterns)...)

end