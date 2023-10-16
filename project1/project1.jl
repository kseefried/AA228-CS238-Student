using Graphs
using Printf

"""
    write_gph(dag::DiGraph, idx2names, filename)

Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
"""
function write_gph(dag::DiGraph, idx2names, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s,%s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end


function compute(infile, outfile)

    # WRITE YOUR CODE HERE
    # FEEL FREE TO CHANGE ANYTHING ANYWHERE IN THE CODE
    # THIS INCLUDES CHANGING THE FUNCTION NAMES, MAKING THE CODE MODULAR, BASICALLY ANYTHING

    function bayesian_score_component(M,alpha)
        # need loggamma function to do this (from SpecialFunctions.jl)
        p = sum(loggamma.(alpha + M))
        p -= sum(loggamma.(alpha))
        p += sum(loggamma.(sum(alpha,dims=2)))
        p-= sum(loggamma.(sum(alpha,dims=2)+sum(M,dims=2)))
    end 

    function bayesian_score(vars, G, D)
        n = length(vars)
        M = statistics(vars, G, D)
        alpha = prior(vars, G)
        return sum(bayesian_score_component(M[i],alpha[i]) for i in 1:n)
    end 

    struct K2search
        ordering::Vector{Int}
    end 

    function fit(method::K2Search, vars, D)
        G = SimpleDiGraph(length(vars))
        for (k,i) in enumerate(method.ordering[2:end])
            y = bayesian_score(vars, G, D)
            for j in method.ordering[1:k]

    end 

end

if length(ARGS) != 2
    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
end

inputfilename = ARGS[1]
outputfilename = ARGS[2]

compute(inputfilename, outputfilename)
