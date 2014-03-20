function [indiv, objEM] = refinement(indiv, data, cfgPrm)
%Refines each individual of the population using EM

if isfield(cfgPrm,'DEBUG')
	DEBUG = cfgPrm.DEBUG;
else
	DEBUG = 0;
end

[nObjects nFeatures] = size(data);

statOpts = statset('MaxIter', cfgPrm.maxEMIter, 'TolFun', 1e-5);

[nClusters means covs mixingCoefficients] = gmm_parameters_from_individual(indiv, nFeatures);

if DEBUG
	fprintf(DEBUG,'#REFINEMENT\nOLD INDIVIDUAL (%d):%s\n',i,info_individual(indiv));
end


%if we find a problem with ill conditioned covariance matrices, restart a new solution
try
	objEM = gmdistribution.fit(data, nClusters, ...
		'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
		'Options', statOpts, 'Regularize', cfgPrm.regV);
catch err
	if DEBUG
		fprintf(DEBUG,'\nProblem with refinement, keeping old solution.\n%s\n',info_individual(indiv))
	end
	objEM = gmdistribution(means, covs, mixingCoefficients);
	return
end

%update individual parameters
indiv.mean(1:nClusters,:) = objEM.mu;
indiv.mixCoef(1:nClusters) = objEM.PComponents;
indiv.fitness = fitnessFunc( cfgPrm.fitnessFName, objEM, nObjects, nClusters, nFeatures );
covs = objEM.Sigma;
for k=1:nClusters
	indiv.covariance(k,:) = squareformSymmetric( covs(:,:,k) );
%	indiv.determinant(k) = det( covs(:,:,k) );
	%storing the squared mahalanobis distance
%	dif = bsxfun(@minus, data, indiv.mean(k,:));
%	indiv.distance(:,k) = sum((dif / covs(:,:,k)) .* dif,2);
end

if DEBUG
	fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(indiv));
end

end


function [nClusters means covs mixingCoefficients] = gmm_parameters_from_individual(indiv,nFeatures)
	nClusters = indiv.nClusters;
	covs = zeros( nFeatures, nFeatures, nClusters );
	means = indiv.mean(1:nClusters,:);

	%create gmdistribution
	for k=1:nClusters
		covs(:,:,k) = squareformSymmetric( indiv.covariance(k,:) );
	end
	mixingCoefficients = indiv.mixCoef(1:nClusters);
end
