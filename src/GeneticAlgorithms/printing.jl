
# ELEMENTS ######################

function Base.show(io::IO, sol::Sol)
    #print(io,"($(sol.z), $(sol.x))")
    #print(io,sol.x)
    print(io,"($(sol.z), $(sol.rank), $(sol.crowding), $(sol.hv_contrib))")
end

# VECTORS ######################

function Base.show(io::IO, vec::Vector{Sol})
    print(io,"[\n")
    for i in 1:length(vec)
        print(io,vec[i])
        print(io,",\n")
    end
    print(io,"]\n")
end
