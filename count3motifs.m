function [total_motifs, single_vertex_motifs] = count3motifs(g)

% Counts 3 motifs in binary directed graph in the scales of single vertex
% and the overall network.
%
% Receives:
%       g - matrix Mx2 - List of M directed links. Each link is weighted
%       equally. In the preprocessing the graph is "squeezed" - if in the
%       input edge list there are uniquely N vertices but some are numbered
%       to values higher than N, the vertices are renumbered to the range
%       of 1...N.
%       
% Returns:
%       single_vertex_motifs - a matrix of N*13 - where N is the number of
%       vertices in g and 13 is the number of different directed motifs of
%       order 4. The i-th row points out the amount of motifs the i-th 
%       vertex participates with for each of the 199 types.
%		total_motifs - 2*13 - where the first column is motif identifier
%		and the second column is number of occurences of this motif overall
%		the network.

addpath(genpath(pwd));
Graph = ObjectCreateGraph(g);
total_motifs = mex3motifs(Graph);
total_motifs = total_motifs(2:end,:);
single_vertex_motifs = load('example3.txt');
single_vertex_motifs = single_vertex_motifs(2:end, 2:end);
delete('example3.txt');
if exist('debugeer.txt','file')
    delete('debugeer.txt');
end


end

