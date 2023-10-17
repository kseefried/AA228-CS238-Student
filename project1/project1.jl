using Graphs
using Printf
using SpecialFunctions
using DataFrames
using LinearAlgebra
using TikzGraphs
using TikzPictures
using CSV

struct Variable
    name::Symbol
    r::Int # number of possible values 
end 

function sub2ind(siz, x)
    k = vcat(1,cumprod(siz[1:end-1]))
    return dot(k, x .-1) + 1
end 

# Extracts statistics and counts from a discrete data set Dict
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
    p = sum(loggamma.(alpha + M))
    p -= sum(loggamma.(alpha))
    p += sum(loggamma.(sum(alpha,dims=2)))
    p-= sum(loggamma.(sum(alpha,dims=2)+sum(M,dims=2)))
end 

# Calculate Total Bayesian Score
function bayesian_score(vars, G, D)
    n = length(vars)
    M = statistics(vars, G, D)
    alpha = prior(vars, G)
    return sum(bayesian_score_component(M[i],alpha[i]) for i in 1:n)
end 

struct K2Search
    ordering::Vector{Int}
end 

function fit(method::K2Search, vars, D)
    G = SimpleDiGraph(length(vars))
    for (k,i) in enumerate(method.ordering[2:end])
        y = bayesian_score(vars, G, D)
        while true 
            y_best, j_best = -Inf, 0
            for j in method.ordering[1:k]
                if !has_edge(G, j, i)
                    add_edge!(G,j,i)
                    y_1 = bayesian_score(vars, G, D)
                    if y_1 > y_best
                        y_best, j_best = y_1, j
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

# Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
function write_gph(dag::DiGraph, idx2names, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s,%s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end

function compute(infile, outfile1)  
        # Get data from CSV file and put it into Variable with a name and max possible number
        data = CSV.File(infile) |> DataFrame
        cols = names(data)
        D = data[2:end,:]
        x = length(names)
        vars = Variable[]

        # Add all variable names and maximum values to vars
        for h in 1:x
            vals = D[:,h]
            max_val = maximum(vals)
            append = Variable(Symbol(names[h]),max_val)
            push!(vars, append)
        end 

        # Run K2 search using extracted data, calculate Bayesian score, write .gph
        fit(method::K2Search,vars,D)
        score = bayesian_score(vars,G,D)
        println(score)
        write_gph(G,idxenames,outfile1)     
end 

inputfilename = "data/small.csv"
outputfilename1 = "small.gph"

@time begin
compute(inputfilename, outputfilename1)
end 

# Using K2 search to find the structure of the data 
# come up with a search method to find the underlying 
# stores the best graph so far after so many iterations, keep track
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
# Function to write the image of parents 