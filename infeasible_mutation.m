function [P] = infeasible_mutation(P, maxClusters, data, constraints)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population

global DEBUG;

if ischar(P) && strcmp(P,'debug')
	unittests();
	return
end


[idxVio,totPenalty,penaltiesByCon]=compute_penalty(P, constraints, data)
[nChunklets,chunklets] = generate_chunklets(conGraph);

for i=1:length(P)

	numClusters = P(i).numClusters;
	penalties = penaltiesByCon(i,:);
	penalties = penalties./sum(penalties);
	maxClustersToCreate = maxClusters - numClusters;
	constraintsToFix = randi(maxClustersToCreate,1) 
	chosen = roulette(penalties,constraintsToFix)

	pdfs = computePosterior(P(i), data);
	for c=1:constraintsToFix
		con = constraints(chosen(c),:)
		[po1,ko1] = pdfs(con(1),:)
		%TODO Terminar aqui
	end
	pElim = (numClusters-2)/(maxClusters-2);
	if (rand() > pElim)
		mutOpToApply = 'create';
	else
		mutOpToApply = 'elim';
	end

	if DEBUG
		fprintf(DEBUG,'\n\n\n-----%s\n\n#MUTATION\nOLD INDIVIDUAL (%d):%s\n',mutOpToApply,i,info_individual(P(i)));
	end


	[posterior gauss] = computePosterior( P(i), data );
	oldIndiv = P(i);
	if strcmp(mutOpToApply,'elim')
		%remove clusters
		lnTerms = bsxfun(@plus, log(eps+P(i).mixCoef(1:numClusters)), log(eps+gauss));
		%transform in positive so that higher values will have more probability to mutate
		valuesMutationClusters = -1*sum(posterior .* lnTerms, 1);
		[valuesSorted idxSorted] = sort( valuesMutationClusters );
		valuesMutationClusters( idxSorted ) = 1:numClusters;

		probs = valuesMutationClusters ./ sum(valuesMutationClusters);
		for c=1:length(unique(P(i).classOfCluster))
			if sum(P(i).classOfCluster == c) == 1
				probs(P(i).classOfCluster == c) = 0
			end
		end
		z = randi([1 numClusters-2]);
		%chosen are the clusters selected for removal
		chosen = roulette_without_reposition( probs, z );
		survivors = setdiff(1:numClusters, chosen);
		P(i).mean = P(i).mean(survivors,:);
		P(i).covariance = P(i).covariance(survivors,:);
		P(i).mixCoef = P(i).mixCoef(survivors);
		%re-normalize
		P(i).mixCoef = P(i).mixCoef ./ nansum(P(i).mixCoef);
		P(i).numClusters = length(survivors);
	else
		%create new clusters
		logsPost = log2(eps+posterior);
		objEntropies = -sum(posterior.*logsPost,2);

		[valuesSorted idxSorted] = sort( objEntropies );
		objEntropies( idxSorted ) = 1:length(objEntropies);

		probs = objEntropies./sum(objEntropies);
		z = randi([1 maxClusters-numClusters]);
		%chosen are the objects that will be used to create clusters
		chosen = roulette_without_reposition( probs, z );
		variances = var(data);
		P(i).numClusters = numClusters + z;
		%recover the clusters that were most probable to generate these points
		[~,oldClusters] = max(posterior(chosen,:), [], 2);
		for nc=1:z
			P(i).mean(numClusters+nc,:) = data(chosen(nc),:);
			P(i).covariance(numClusters+nc,:) = squareformSymmetric( 0.1*diag(variances) );
			P(i).mixCoef(numClusters+nc) = P(i).mixCoef(oldClusters(nc))/2;
			P(i).mixCoef(oldClusters(nc)) = P(i).mixCoef(oldClusters(nc))/2;
			P(i).classOfCluster(numClusters+nc,:) = pickClass(chosen(nc), data, chunklets)
		end
	end

	if DEBUG
		assert(P(i).numClusters >= 2, 'Too few clusters')
		assert(P(i).numClusters <= maxClusters, 'Too many clusters')
		fprintf(DEBUG,'#MUTATION\nNEW INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
	end


end

end

function [closestClass] = pickClass(idxObj, data, chunklets)
	if chunklets(idxObj) > 0
		closestClass = chunklets(idxObj);
	else
		labeled = data(chunklets>0,:);
		classes = chunklets(chunklets>0);
		closestClass = NaN;
		dists = pdist2(data(idxObj,:), labeled, 'euclidean');
		dists = max(dists)-dists;
		dists = dists/sum(dists);
		chosen = roulette(dists, 1);
		closestClass = classes(chosen);
	end
end

function unittests
	testPickClass;
end

function testPickClass
	chunklets = [1 0 1 0 0 2 2];
	data = [ 1 2; 2 2; 1 2; 4 4; 80 80; 100 90; 120 100];
	assert(pickClass(1, data, chunklets) == 1)
	assert(pickClass(2, data, chunklets) == 1)
	assert(pickClass(3, data, chunklets) == 1)
	assert(pickClass(4, data, chunklets) == 1)
	assert(pickClass(5, data, chunklets) == 2)
	assert(pickClass(6, data, chunklets) == 2)
	assert(pickClass(7, data, chunklets) == 2)
end

