load data/cctoolbox.mat
addpath('libs/CCToolbox');
setcctpath;
M = curve_clust(trajs,Options);
showmodel(M,trajs)
