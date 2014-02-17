function [P] = refinement(P, data, cfgPrm)
%Refines each individual of the population using EM

global EMSteps;
global DEBUG;

[nObjects nFeatures] = size(data);

statOpts = statset('MaxIter', cfgPrm.maxEMIter, 'TolFun', 1e-5);
%regularize value
regV = 1e-5;
fitnessFName='mdl';
for i=1:length(P)


	[nClusters means covs mixingCoefficients] = gmm_parameters_from_individual(P(i));

	if DEBUG
		fprintf(DEBUG,'#REFINEMENT\nOLD INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
	end


	%if we find a problem with ill conditioned covariance matrices, restart a new solution
	try
		objEM = gmdistribution.fit(data, nClusters, ...
			'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
			'Options', statOpts, 'Regularize', regV);
	catch err
		if DEBUG
			fprintf(DEBUG,'\nProblem with refinement, keeping old solution.\n%s\n',info_individual(P(i)))
		end
		continue
	end

	EMSteps = EMSteps + objEM.Iters;

	%update individual parameters
	P(i).mean(1:nClusters,:) = objEM.mu;
	P(i).mixCoef(1:nClusters) = objEM.PComponents;
	P(i).fitness = fitnessFunc( fitnessFName, objEM, nObjects, nClusters, nFeatures );
	covs = objEM.Sigma;
	for k=1:nClusters
		P(i).covariance(k,:) = squareformSymmetric( covs(:,:,k) );
		P(i).determinant(k) = det( covs(:,:,k) );
		%storing the squared mahalanobis distance
		dif = bsxfun(@minus, data, P(i).mean(k,:));
		P(i).distance(:,k) = sum((dif / covs(:,:,k)) .* dif,2);
	end

	if DEBUG
		fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(P(i)));
	end

end

end

function [nClusters means covs mixingCoefficients] = gmm_parameters_from_individual(indiv)
	nClusters = indiv.nClusters;

	covs = zeros( nFeatures, nFeatures, nClusters );
	means = indiv.mean(1:nClusters,:);

	%create gmdistribution
	for k=1:nClusters
		%adding regularization value to avoid ill conditioning
		covs(:,:,k) = squareformSymmetric( P( i ).covariance(k,:) ) + eye(nFeatures)*regV;
	end
	mixingCoefficients = indiv.mixCoef(1:nClusters);
end
