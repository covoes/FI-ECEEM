function [Pfeas Pinfeas] = initialize_with_constraints(data, conGraph, maxClusters, 
	sizePopulationFeasible, sizePopulationInfeasible, maxKMSIter)
%Initialize individuals validating them with the constraints
	
	Pfeas = []
	Pinfeas = []

	[nChunklets,chunklets] = graphconncomp(conGraph) 
	[numObjects numFeatures] = size(data);

	empty_individual = create_empty_individual(nObjects, nFeatures, maxClusters)

	Pfeas = repmat(empty_individual, [1 sizePopulationFeasible]);
	Pinfeas = repmat(empty_individual, [1 sizePopulationInfeasible]);

	opts = statset('MaxIter',maxKMSIter);

	if sizePopulationFeasible > 1
		numClusters = floor(linspace(2,maxClusters,sizePopulationFeasible));
	else
		numClusters = randi([2 maxClusters],1);
	end

	while len(Pfeas) < sizePopulationFeasible 		
		
		[idx, clusters] = kmeans(data, numClusters(i), 'Start', init, 'EmptyAction', 'singleton',...
			                       'Options', opts ); 
		P(i).numClusters = numClusters(i);
		for k=1:numClusters(i)
			P(i).mean(k,:) = clusters(k,:);
			idtmp = find(idx == k);
			%the duplicate is because of the case nObjInCluster=1
			covariance = cov([data(idtmp,:); data(idtmp,:)]);
			P(i).covariance(k,:) = squareformSymmetric( covariance );
			P(i).mixCoef(k) = length(idtmp)/numObjects;
		end



for i=1:sizePopulation
		%Tries to generate more solutions than needed in most cases, this may be faster than trying 
		%one by one due to the validity of solutions
		P = initialize(data, maxClusters, sizePopulationFeasible + sizePopulationInfeasible, maxKMSiter)

		for i=1:length(P)
			if is_valid_solution(P(i), constraints)
				Pfeas(length(Pfeas)+1) = P(i)
			else
				PInfeas(length(Pinfeas)+1) = P(i)
		end
or len(Pinfeas) < sizePopulationInfeasible

end

function individual = create_empty_individual(nObjects, nFeatures, maxClusters)
		individual = struct( 'mean', NaN([maxClusters numFeatures]), ...
					 'covariance', repmat(squareformSymmetric(NaN(numFeatures))', maxClusters, 1), ...
					 'mixCoef', NaN([1 maxClusters]), ...
					 'distance', NaN([numObjects maxClusters]), ...
					 'determinant', NaN([1 maxClusters]), ...
					 'numClusters', NaN, ...
					 'fitness', NaN);

