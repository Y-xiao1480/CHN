function [centroid, result] = IsoData(data,method, k , it, minimum_n, maximum_variance, minimum_d)
addpath('.\lib');
addpath('.\tool');


if(strcmp(method,'ISODATA') || strcmp(method,'isodata'))
    desired_k=k;  % desired number of classes
    iteration=it;  % maximum iteration time
    minimum_n;  % minimum number of samples in one class
    maximum_variance;  % maximum allowed variance of samples in one class
    minimum_d;  %  minimum distance between two classes
    [centroid, result] = ISODATA(data, iteration, desired_k, minimum_n, maximum_variance, minimum_d);
end