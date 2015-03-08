function [Graph,varargout] = GraphRemoveReciprocalLinks(Graph)
% removes reciprocal links. Of the two reciprocal links, the one with i<j is left. Possible loops (i=j) are preserved
%
% Receives:
%   Graph - structure - Graph. See GraphLoad.
% 
% Returns:
%   Graph - structure - Graph. See GraphLoad.
%   varargout
%       - (optional) - the indecis (in the original graph) of the removed links. 
%
% See Also:
%   GraphSortLinks,GraphSortLinks
%{
    Graph = GraphLoadSample('poisson','N',100,'p',0.05);
    [Graph , removed_links_indices]= GraphRemoveReciprocalLinks(Graph);

%}

narginchk(1,1);
nargoutchk(0,2);

MaxNodeID = max(max(Graph.Data(:,1:2)))+1; 
[~, originally_undirected_ai] = setdiff( Graph.Data(:,1)+ MaxNodeID*Graph.Data(:,2), Graph.Data(:,2)+ MaxNodeID*Graph.Data(:,1));
[~, reciprocal_ai, ~] = intersect( Graph.Data(:,1)+ MaxNodeID*Graph.Data(:,2), Graph.Data(:,2)+ MaxNodeID*Graph.Data(:,1));

ReciprocalLinks = Graph.Data(reciprocal_ai,:); 
if nargout>1
    varargout{1}=reciprocal_ai( ReciprocalLinks(:,1)>ReciprocalLinks(:,2) ); 
end
Graph.Data = [ Graph.Data(originally_undirected_ai,:); ReciprocalLinks(ReciprocalLinks(:,1)<=ReciprocalLinks(:,2),:)];
Graph = GraphSortLinks(Graph);