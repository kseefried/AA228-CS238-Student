using Graphs
using Printf
using SpecialFunctions
using DataFrames
using LinearAlgebra
using TikzGraphs
using TikzPictures
using CSV

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

struct Variable
    name::Symbol
    r::Int # number of possible values 
end 

function extract_data(filename)
    # takes in a file name 
    # extracts the variable names and the data structure 
    data = CSV.File(filename) |> DataFrame
    df = DataFrame(data)
    vars = df[1,:]
    D = df[2:end,:]
    return vars
    return D
end 

function sub2ind(siz, x)
    k = vcat(1,cumprod(siz[1:end-1]))
    return dot(k, x .-1) + 1
end 

# extracts statistics or counts from a discrete data set Dict
# Define all functions before, then run functions that calls functions
function statistics(vars, G,D::Matrix{Int})
    n = size(D,1) 
    r = [vars[i].r for i in 1:n]
    q = [prod([r[j] for j in inneighbors(G,i)]) for i in 1:n]
    M = [zeros(q[i], r[i]) for i in 1:n]
    for o in eachcol(D)
        for i in 1:n
            k = o[i]
            parents = inneighbors(G,i)
            j = 1
            if !is empty(parents)
                j = sub2ind(r[parents],o[parents])
            end 
            M[i][j,k]+= 1.0
        end 
    end 
    return M
end 

# Generate Uniform Prior
function prior(vars, G) 
    n = length(vars)
    r = [vars[i].r for i in 1:n]
    q = [prod([r[j] for j in inneighbors(G,i)]) for i in 1:n]
    return [ones(q[i],r[i]) for i in 1:n]
end 

# Calculate Bayesian Score Component
function bayesian_score_component(M,alpha)
    # need loggamma function to do this (from SpecialFunctions.jl)
    p = sum(loggamma.(alpha + M))
    p -= sum(loggamma.(alpha))
    p += sum(loggamma.(sum(alpha,dims=2)))
    p-= sum(loggamma.(sum(alpha,dims=2)+sum(M,dims=2)))
end 

# Calculate total bayesian score
function bayesian_score(vars, G, D)
    n = length(vars)
    M = statistics(vars, G, D)
    alpha = prior(vars, G)
    return sum(bayesian_score_component(M[i],alpha[i]) for i in 1:n)
end 

# Using K2 search to find the structure of the data 
# come up with a search method to find the underlying 
# stores the best graph so far after so many iterations, keep track

struct K2search
    ordering::Vector{Int}
end 

function fit(method::K2Search, vars, D)
    G = SimpleDiGraph(length(vars))
    for (k,i) in enumerate(method.ordering[2:end])
        y = bayesian_score(vars, G, D)
        while true 
            for j in method.ordering[1:k]
                if !has_edge(G, j, i)
                    y' = bayseian_score(vars, G, D)
                    if y' > y_best
                        y_best, j_best = y', j
                    end 
                    rem_edge!(G, j, i) 
                end 
            end 
            if y_best > y
                y = y_best 
                add_edge!(G, j_best, i)
            else 
                break
            end 
        end 
    end 
    return G
end 


function compute(infile, outfile)
    extract_data(infile)
# define variable struct
# need csv file variables 
# need values taken on my the values 
# get counts for everything, extracting stats from everything 
# define the prior, assume uniform Prior
# method for Bayesian score 
# more creative -- different types of searches
# -- tweak things, inject randomness
# check bayesian score calculation using the example 

# how to print out graph

# need bayesian score, .gph, and .pdf 
    
end

if length(ARGS) != 2
    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
end

inputfilename = ARGS[1]
outputfilename = ARGS[2]


# Function to write the image of parents 

compute(inputfilename, outputfilename)
