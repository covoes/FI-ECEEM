function [indiv] = feasible_mutation(indiv, gmmObj, sharedData, configPrm)
%FEASIBLE_MUTATION Performs mutation on the indivuals in the population

	if ischar(indiv) && strcmp(indiv,'debug')
		unittests()
		return
	end

	data = sharedData.data;
	chunklets = sharedData.chunklets;

	nClusters = indiv.nClusters;

	pElim = (nClusters-2)/(configPrm.maxClusters-2);
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
		[valuesSorted idxSorted] = sort( valuesMutationClusters );
		valuesMutationClusters( idxSorted ) = 1:nClusters;
		probs = valuesMutationClusters ./ sum(valuesMutationClusters);
		if all(probs==0)
			return
		end	
		z = randi([1 nClusters-2]);
		%chosen are the clusters selected for removal
		chosen = roulette_without_reposition( probs, z );
		chosen = filter_singleton_classes_clusters(chosen, indiv.classOfCluster);
		if isequal(chosen,[])
			return
		end
		survivors = setdiff(1:nClusters, chosen);
		indiv.mean = indiv.mean(survivors,:);
		indiv.covariance = indiv.covariance(survivors,:);
		indiv.mixCoef = indiv.mixCoef(survivors);
		indiv.classOfCluster = indiv.classOfCluster(survivors);
		%re-normalize
		indiv.mixCoef = indiv.mixCoef ./ nansum(indiv.mixCoef);
		indiv.nClusters = length(survivors);
	else
		%create new clusters
		logsPost = log2(eps+posterior);
		objEntropies = -sum(posterior.*logsPost,2);

		[valuesSorted idxSorted] = sort( objEntropies );
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


end

function [filtered] = filter_singleton_classes_clusters(chosen, classOfCluster)
	filtered = chosen;
	classesOfChosen = classOfCluster(chosen);
	for c=1:max(classOfCluster)
		nClustersOfClass = sum(classOfCluster == c);
		chosenClustersOfClass = find(classesOfChosen == c);
		if nClustersOfClass == 1
			filtered(chosenClustersOfClass) = -1;
		elseif nClustersOfClass == length(chosenClustersOfClass)
			filtered(chosenClustersOfClass) = -1;
			filtered(length(filtered)+1) = randsample(chosen(chosenClustersOfClass),1);
		end
	end
	filtered = filtered(filtered>-1);
end

function unittests
	test_filterSingletonClassesClusters
end

function test_filterSingletonClassesClusters
	chosen = [ 2 4 ];
	classOfClusters = [ 1 2 1 2];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, 2) || isequal(ret,4));
	
	chosen = 2 ;
	classOfClusters = [ 1 2 1 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, []))

	chosen = [ 1 4 ];
	classOfClusters = [ 1 2 1 3 4 3 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue(isequal(ret, [1 4]))

	chosen = [ 1 3 4 6 ];
	classOfClusters = [ 1 2 1 3 4 3 ];
	ret = filter_singleton_classes_clusters(chosen, classOfClusters);
	assertTrue((ret(1) == 1 || ret(1) == 3) && (ret(2) == 4 || ret(2) == 6))
  assertTrue(length(ret) == 2)	
end


