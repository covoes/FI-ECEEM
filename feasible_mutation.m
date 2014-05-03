function [indiv,gmmObj] = feasible_mutation(indiv, gmmObj, sharedData, configPrm)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population
	data = sharedData.data;
	chunklets = sharedData.chunklets;

	nClusters = indiv.nClusters;

	pElim = (nClusters-configPrm.minClusters)/(configPrm.maxClusters-configPrm.minClusters);
	if (rand() > pElim)
		mutOpToApply = 'create';
	else
		mutOpToApply = 'elim';
	end

	posterior = gmmObj.posterior;
	gauss = gmmObj.pdf;
	oldIndiv = indiv;
	if strcmp(mutOpToApply,'elim')
		%remove clusters
		lnTerms = bsxfun(@plus, log(eps+indiv.mixCoef(1:nClusters)), log(eps+gauss));
		%transform in positive so that higher values will have more probability to mutate
		valuesMutationClusters = -1*sum(posterior .* lnTerms, 1);
		[~, idxSorted] = sort( valuesMutationClusters );
		valuesMutationClusters( idxSorted ) = 1:nClusters;
		probs = valuesMutationClusters ./ sum(valuesMutationClusters);
		if all(probs==0)
			return
		end	
		z = randi([1 nClusters-configPrm.minClusters]);
		%chosen are the clusters selected for removal
		chosen = roulette_without_reposition( probs, z );
		indiv = remove_clusters(chosen, indiv);
	else
		%create new clusters
		logsPost = log2(eps+posterior);
		objEntropies = -sum(posterior.*logsPost,2);

		[~, idxSorted] = sort( objEntropies );
		objEntropies( idxSorted ) = 1:length(objEntropies);

		probs = objEntropies./sum(objEntropies);
		z = randi([1 configPrm.maxClusters-nClusters]);
		%chosen are the objects that will be used to create clusters
		chosen = roulette_without_reposition( probs, z );
		%recover the clusters that were most probable to generate these points
		oldClusters = gmmObj.clusterLabels;
		for nc=1:z
			indiv = create_cluster(indiv, chosen(nc), data, chunklets, oldClusters);
		end
	end


	if configPrm.DEBUG
		fprintf(configPrm.DEBUG,'\n\n\n-----%s\n\n#MUTATION\nOLD INDIVIDUAL :%s\n',...
			mutOpToApply,info_individual(oldIndiv));
		fprintf(configPrm.DEBUG,'#MUTATION\nNEW INDIVIDUAL :%s\n',info_individual(indiv));
		assert(indiv.nClusters >= 2, 'Too few clusters')
		assert(indiv.nClusters <= configPrm.maxClusters, 'Too many clusters')
		assert(isequal(sort(unique(indiv.classOfCluster(:)))', 1:max(chunklets)))
	end

	if indiv.nClusters > 2
		indiv = clean_solution(indiv,sharedData, configPrm);
	end
	[indiv,gmmObj] = refinement(indiv, sharedData, configPrm);

end

