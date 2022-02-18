function [centroid, result] = kmeans_pro(data, method, k,it)
addpath('.\lib');
addpath('.\tool');


% Kmeans++
if(strcmp(method,'kmeans++') || strcmp(method,'kmeanspp')) 
    k
    iteration = it;
    [centroid, result] = Kmeanspp(data, k, iteration);
end