function [Pfeas Pinfeas] = initialize_with_constraints(sharedData, configPrm)
%Initialize individuals validating them with the constraints

  if ischar(sharedData) && strcmp(sharedData,'debug')
		unittests()
		return
	end

	Pfeas = [];
	Pinfeas = [];
	sizePopFeasible = configPrm.sizePopulationFeasible;
	sizePopInfeasible = configPrm.sizePopulationInfeasible;
	nClusters = generate_nclusters(sizePopFeasible,configPrm);

	data = sharedData.data;
	nChunklets = sharedData.nChunklets;
	chunklets = sharedData.chunklets;

	nTry = 0;
	while length(Pfeas) < sizePopFeasible && configPrm.maxInitTries > nTry
		nTry = nTry + 1;
		k = nClusters(mod(nTry,sizePopFeasible)+1);
		classOfCluster = spread_chunklets_among_clusters(nChunklets, k);
		seeds = sample_initial_seeds_off_chunklets(chunklets, classOfCluster);
		initMatrix = data(seeds,:);
		[individual,gmmObj] = initialize_gmm(sharedData, k, initMatrix, configPrm, classOfCluster);
		assert (isequal(sort(unique(individual.classOfCluster(:)))', 1:max(chunklets)))
		[Pfeas Pinfeas] = insert_individual_correct_pool(individual,gmmObj, Pfeas, Pinfeas);
	end

	nTry = 0;
	nClasses = max(chunklets);
	nClusters = generate_nclusters(sizePopInfeasible,configPrm);
	while length(Pinfeas) < sizePopInfeasible && configPrm.maxInitTries > nTry
		nTry = nTry + 1;
		k = nClusters(mod(nTry,sizePopInfeasible)+1);
		initMatrix = data(randsample(size(data,1),k),:);
		classOfCluster = randsample(nClasses, nClasses);
		if k > nClasses
			classOfCluster(nClasses+1:k) = randi(nClasses, [1 k-nClasses]);
		end	
		[individual,gmmObj] = initialize_gmm(sharedData, k, initMatrix, configPrm, classOfCluster);

		assert (isequal(sort(unique(individual.classOfCluster(:)))', 1:max(chunklets)))
		[Pfeas Pinfeas] = insert_individual_correct_pool(individual,gmmObj, Pfeas, Pinfeas);
	end

	Pfeas = Pfeas(1:min(length(Pfeas),sizePopFeasible));
	Pinfeas = Pinfeas(1:min(length(Pinfeas),sizePopInfeasible));
end

function [individual,gmmObj] = initialize_gmm(sharedData, k, init,configPrm, classOfCluster)
	opts = statset('MaxIter',configPrm.maxKMSIter);
	[idx, clusters] = kmeans(sharedData.data, k, 'EmptyAction', 'singleton',...
	                       	'Start', init, 'Options', opts);
	individual = gmmFromKmeans(idx, clusters, sharedData.data);
	individual.classOfCluster = classOfCluster;
	[individual,gmmObj] = one_em_step(individual, sharedData, configPrm);
end

function [individual,gmmObj] = one_em_step(individual, sharedData, configPrm)
	[individual,gmmObj] = refinement(individual, sharedData, configPrm);
end

function seeds = sample_initial_seeds_off_chunklets(chunklets, classOfCluster)
	seeds = NaN([1 length(classOfCluster)]);
	for k=1:length(classOfCluster)
		candidates = find(chunklets==classOfCluster(k));
		if length(candidates) == 1
			seeds(k) = candidates;
		else
			seeds(k) = randsample(candidates, 1);
		end
	end
end


function gmm = gmmFromKmeans(idx, clusters, data)
	gmm = create_empty_individual;
	gmm.nClusters = size(clusters,1);
	[nObjects nFeatures] = size(data);
	for k=1:gmm.nClusters
		gmm.mean(k,:) = clusters(k,:);
		idtmp = find(idx == k);
		if length(idtmp) == 1
			covariance = eye(nFeatures)*eps;
		else
			covariance = cov(data(idtmp,:));
		end
		gmm.covariance(k,:) = squareformSymmetric( covariance );
		gmm.mixCoef(k) = length(idtmp)/nObjects;
	end
end


function [classOfCluster] = spread_chunklets_among_clusters(nChunklets, k)
	clustersToMap = k;
	classOfCluster = [];
	while clustersToMap > 0
		nToSample = min([nChunklets, clustersToMap]);
		classOfCluster = [classOfCluster;randsample(nChunklets,nToSample)];
		clustersToMap = clustersToMap - nToSample;
	end
end

function individual = create_empty_individual
		individual = struct( 'mean', [], ...
					 'covariance', [], ...
					 'mixCoef', [], ...
			     'classOfCluster', [], ...
					 'nClusters',NaN, ...
					 'fitness', NaN);
end

function nClusters = generate_nclusters(nIndividuals,cfg)
	if nIndividuals > 1
		nClusters = floor(linspace(cfg.minClusters,cfg.maxClusters,nIndividuals));
	else
		nClusters = randi([cfg.minClusters cfg.maxClusters],1);
	end
end

function unittests
	testSpreadChunklets
	testGenerateNClusters
	testGmmFromKmeans
	testSampleInitialSeedsOfChunklets
	testInitializeWithConstraints
end

function testInitializeWithConstraints
	data = [ 1.1 2; 1.2 3; 2 2; 5.1 5.2; 9 10; 10 11; 21.1 22.4];
	con = [1 3 1;3 5 -1; 4 6 -1; 5 6 1];
	conHard = generate_constraint_graph(con, size(data,1));
	configPrm = struct('minClusters',2,'maxClusters',4, 'maxKMSIter',3, ...
		                 'sizePopulationFeasible', 4, 'sizePopulationInfeasible',4,...
	                   'maxInitTries', 10, 'DEBUG',0,'regV', 1e-2, ...
		                 'maxEMIter', 2, 'fitnessFName', 'mdl');

	rng(42);
	[nChunklets chunklets] = generate_chunklets(conHard);
	shared = struct('data', data, 'constraints', con,'chunklets', chunklets,...
	               	'nChunklets', nChunklets, 'conGraph', conHard);
	[outFeas outInfeas] = initialize_with_constraints(shared, configPrm);
	assertEqual(length(outFeas), configPrm.sizePopulationFeasible)
	assertEqual(length(outInfeas), configPrm.sizePopulationInfeasible)
	assertTrue(all([outFeas(:).nClusters] >= configPrm.minClusters),'Numero de clusters invalido')
  assertTrue(all([outFeas(:).nClusters] <= configPrm.maxClusters),'Numero de clusters invalido')
	assertTrue(all([outInfeas(:).nClusters] >= configPrm.minClusters),'Numero de clusters invalido')
  assertTrue(all([outInfeas(:).nClusters] <= configPrm.maxClusters),'Numero de clusters invalido')
	%TODO reimplementar teste sem necessitar do distance dentro do individuo
%	for i=1:length(outFeas)
%		for c=1:size(con,1)
%			[~,k1] = min(outFeas(i).distance(con(c,1),:));
%			[~,k2] = min(outFeas(i).distance(con(c,2),:));
%			if con(c,3) == 1
%			 	assertEqual(outFeas(i).classOfCluster(k1), outFeas(i).classOfCluster(k2),...
%					'Restrição ML violada em solucão factível')
%			else
%				assertTrue(outFeas(i).classOfCluster(k1) ~= outFeas(i).classOfCluster(k2),...
%					'Restrição CL violada em solucão factível')
%			end
%		end
%	end
%
%	for i=1:length(outInfeas)
%		nvio = 0;
%		for c=1:size(con,1)
%			[~,k1] = min(outFeas(i).distance(con(c,1),:));
%			[~,k2] = min(outFeas(i).distance(con(c,2),:));
%			if con(c,3) == 1 && outFeas(i).classOfCluster(k1) ~= outFeas(i).classOfCluster(k2)
%					nvio = nvio + 1;
%			elseif con(c,3) == -1 && outFeas(i).classOfCluster(k1) == outFeas(i).classOfCluster(k2)
%					nvio = nvio + 1;
%			end
%			assertTrue(nvio==0, 'Solução sem violações no pool de não factíveis')
%		end
%	end
end

function testSpreadChunklets
	out = spread_chunklets_among_clusters(10, 10);
	assertEqual(sort(out)', 1:10);

	out = spread_chunklets_among_clusters(5, 10);
	assertEqual(sort(out(1:5))', 1:5);
	assertEqual(sort(out(6:10))', 1:5);
end

function testGenerateNClusters
	out = generate_nclusters(5,struct('minClusters',2,'maxClusters', 6));
	assertEqual(out, 2:6)
	out = generate_nclusters(3,struct('minClusters',4,'maxClusters', 6));
	assertEqual(out, 4:6)
	out = generate_nclusters(1,struct('minClusters',2,'maxClusters', 10));
	assertEqual(length(out),1)
	assertTrue(all(out >= 2) && all(out <= 10))
end

function testGmmFromKmeans
	clusters =  [1 3; 5 7];
	out = gmmFromKmeans([1 1 1 2 2],clusters, [ 1 2; 1 3; 1 4; 4 6; 6 8]);
	assertEqual(out.nClusters, 2)
	assertEqual(out.mean, clusters)
	assertEqual(out.covariance, [0 0 1; 2 2 2])
	assertEqual(out.mixCoef, [0.6 0.4])

	out = gmmFromKmeans([1 1 1 2],clusters, [ 1 2; 1 3; 1 4; 5 7]);
	assertEqual(out.nClusters, 2)
	assertEqual(out.mean, clusters)
	assertElementsAlmostEqual(out.covariance, [0 0 1; 0 0 0])
	assertElementsAlmostEqual(out.mixCoef, [0.75 0.25])
end


function testSampleInitialSeedsOfChunklets
	out = sample_initial_seeds_off_chunklets([ 1 1 2 2 3], [3 2 1]);
	assertTrue(length(out) == 3)
	assertTrue(out(1) == 5)
	assertTrue(out(2) == 3 || out(2) == 4)
	assertTrue(out(3) == 1 || out(3) == 2)
end
