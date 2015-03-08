function Graph = GraphSortLinks(Graph)
% sorts links in ascending order, by i, j. 
%
% Receives:
%   Graph - structure - Graph. See GraphLoad.
% 
% Returns:
%   Graph - structure - Graph. See GraphLoad.
%   varargout
%       - (optional) - the indecis (in the original graph) of the removed links. 
%
%{
    Graph = GraphLoadSample('poisson','N',100,'p',0.05);
    [Graph , removed_links_indices]= GraphRemoveReciprocalLinks(Graph);
%}


MaxNodeID = max(max(Graph.Data(:,1:2)))+1; 
[~, SO] =  sort(Graph.Data(:,1)*MaxNodeID+Graph.Data(:,2));
Graph.Data = Graph.Data(SO,:);