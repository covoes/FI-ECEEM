function s=info_individual(ind)
%INFO_INDIVIDUAL Return a string with all the information on the current individual.
%
%Used for debugging purposes

%print just the 5 first object distances.. 
s = sprintf('\n\tMEAN: %s\n\tCOVARIANCE %s\n\tMIXCOEF %s\n\tNUMCLUSTERS %d\n\tFITNESS %.4f\n\tCLASSES %s\n',mat2str(ind.mean,4), mat2str(ind.covariance,4), mat2str(ind.mixCoef,4), ind.nClusters, ind.fitness, mat2str(ind.classOfCluster,0));

