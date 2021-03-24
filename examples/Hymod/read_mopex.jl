using Dates

function read_mopex(;url=nothing)
    if !startswith(url, "file://")
        error("only local file:// urls are supported")
    end
    _, path = split(url, "file://")
    nlines = 0 
    open(path, read=true) do fh
        # count the number of datapoints
    	nlines = countlines(fh)
	seekstart(fh)

	# preallocate the data vectors
	# given the number of datapoints in series
	datetime = Vector{DateTime}(undef, nlines)
	precip = Vector{Float64}(undef, nlines)
	pet = Vector{Float64}(undef, nlines)
	flow = Vector{Float64}(undef, nlines)
	temp = Vector{Float64}(undef, nlines)
		
	# replace float values with -99.0 set value to NaN
	replace_nan(v::Float64) = v == -99.0 ? NaN : v

	# parse each line
	for (nline, line) in enumerate(readlines(fh))
	    # parse date
	    year = parse(Int, line[1:4])  
	    month = parse(Int, line[5:6])
	    day = parse(Int, line[7:8])
	    datetime[nline] = DateTime(year, month, day)
	    # parse numeric values
	    _precip, _pet, _flow, _max_temp, _min_temp = 
	    	map(x -> parse(Float64, x), split(line[9:end]))
	    # precip mm / day
	    precip[nline] = replace_nan(_precip)
	    # pet mm / day
	    pet[nline] = replace_nan(_pet)
	    # flow mm 
	    flow[nline] = replace_nan(_flow)
	    # mean daily temp C
	    temp[nline] = (replace_nan(_min_temp) + replace_nan(_max_temp)) / 2.0
    	end	
    	return (datetime=datetime, precip=precip, pet=pet, flow=flow, temp=temp)
    end
end
