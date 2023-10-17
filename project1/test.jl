using Graphs  # for DiGraph and add_edge!
using TikzGraphs   # for TikZ plot output
using TikzPictures # to save TikZ as PDF

g = DiGraph(2) # create a directed graph
add_edge!(g, 1, 2) # add edge from node 1 to node 2

p = plot(g, ["First", "Second"]) # create TikZ plot with labels
save(PDF("graph.pdf"), p) # save TikZ as PDF