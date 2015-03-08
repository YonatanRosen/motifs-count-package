load c_elegans;
addpath(genpath(pwd));
Graph = ObjectCreateGraph(g);
[total_motifs3, single_vertex_motifs3] = count3motifs(g);
[total_motifs4, single_vertex_motifs4] = count4motifs(g);



