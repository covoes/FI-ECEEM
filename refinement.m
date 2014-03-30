function [indiv, gmmObj] = refinement(indiv, sharedData, cfgPrm, doNotRunEM)
%Refines each individual of the population using EM
%
% If doNotRunEM equals 1, only construct the gmmObj, no optimization is performed.

if isfield(cfgPrm,'DEBUG')
	DEBUG = cfgPrm.DEBUG;
else
	DEBUG = 0;
end

if nargin < 4
	doNotRunEM = 0;
end

data = sharedData.data;
constraints = sharedData.constraints;

[nObjects nFeatures] = size(data);

[nClusters means covs mixingCoefficients objEM] = gmm_parameters_from_individual(indiv, nFeatures);

if doNotRunEM
	gmmObj = gera_struct_gmm_obj() ;
	return
end

statOpts = statset('MaxIter', cfgPrm.maxEMIter, 'TolFun', 1e-5);

if DEBUG
	fprintf(DEBUG,'#REFINEMENT\nOLD INDIVIDUAL :%s\n',info_individual(indiv));
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
	gmmObj = gera_struct_gmm_obj();
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

gmmObj = gera_struct_gmm_obj() ;


function [gmmObj] = gera_struct_gmm_obj
	[idx,nlogl,post,logpdf,mahalad] = cluster(objEM, data);
	pdf = zeros(size(data,1), objEM.NComponents);
	for kk=1:objEM.NComponents
		pdf(:,kk) = mvnpdf(data,objEM.mu(kk,:),objEM.Sigma(:,:,kk));
	end
  [isInfeasible,totPenalty,penaltyByCon] = compute_penalty(indiv, data, constraints, pdf);
	gmmObj = struct('modelo', objEM, 'posterior', post, 'pdf', pdf, ...
		              'clusterLabels', idx, 'classLabels', indiv.classOfCluster(idx), ...
		              'isFeasible', ~isInfeasible, 'totPenalty', totPenalty, ...
		              'penalties', penaltyByCon);
end

end


function [nClusters means covs mixingCoefficients,objEM] = gmm_parameters_from_individual(indiv,nFeatures)
	nClusters = indiv.nClusters;
	covs = zeros( nFeatures, nFeatures, nClusters );
	means = indiv.mean(1:nClusters,:);

	%create gmdistribution
	for k=1:nClusters
		covs(:,:,k) = squareformSymmetric( indiv.covariance(k,:) );
	end
	mixingCoefficients = indiv.mixCoef(1:nClusters);
	objEM = gmdistribution(means,covs,mixingCoefficients);
end
