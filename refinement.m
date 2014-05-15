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
chunklets = sharedData.chunklets;

[nObjects nFeatures] = size(data);

[nClusters means covs mixingCoefficients objEM] = ...
           gmm_parameters_from_individual(indiv, nFeatures, cfgPrm.regV);

gmmObj = [];
if doNotRunEM
	update_parameters();
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
		fprintf(DEBUG,'\nProblem with refinement, keeping old solution.\n%s\n',info_individual(indiv));
		fprintf(DEBUG,err);
	end
%	gmmObj = update_parameters();
%	return
end


update_parameters() ;

if DEBUG
	fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(indiv));
end



function update_parameters
	%update individual parameters
	indiv.mean(1:nClusters,:) = objEM.mu;
	indiv.mixCoef(1:nClusters) = objEM.PComponents;
	covs = objEM.Sigma;
	pdf = zeros(size(data,1), nClusters);
	for k=1:nClusters
		indiv.covariance(k,:) = squareformSymmetric( covs(:,:,k) );
		%[T,P]=cholcov(objEM.Sigma(:,:,kk),0)
		pdf(:,k) = mvnpdf(data,indiv.mean(k,:),...
			                 covs(:,:,k)+eye(size(data,2))*cfgPrm.regV);
	end

	[idx,nlogl,post] = cluster(objEM, data);
  [isInfeasible,totPenalty,penaltyByCon] = compute_penalty(indiv, chunklets, idx, post);
	gmmObj = struct('modelo', objEM, 'posterior', post, 'pdf', pdf, ...
		              'clusterLabels', idx, 'classLabels', indiv.classOfCluster(idx), ...
		              'isFeasible', ~isInfeasible, 'penalties', penaltyByCon);
	indiv.totPenalty = totPenalty;
	indiv.fitness = NaN;
	if ~isInfeasible
		indiv.fitness = fitnessFunc( cfgPrm.fitnessFName, nlogl, nObjects, nClusters, nFeatures );
	end
end

end


function [nClusters means covs mixingCoefficients,objEM] = ...
	                      gmm_parameters_from_individual(indiv,nFeatures,regV)
	nClusters = indiv.nClusters;
	covs = zeros( nFeatures, nFeatures, nClusters );
	means = indiv.mean(1:nClusters,:);

	%create gmdistribution
	for k=1:nClusters
		covs(:,:,k) = squareformSymmetric( indiv.covariance(k,:) ) + eye(nFeatures)*regV;
	end
	mixingCoefficients = indiv.mixCoef(1:nClusters);
	objEM = gmdistribution(means,covs,mixingCoefficients);
end
