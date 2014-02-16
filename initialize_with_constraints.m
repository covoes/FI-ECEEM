function [Pfeas Pinfeas] = initialize_with_constraints(data, conGraph, configPrm)
%Initialize individuals validating them with the constraints

  if ischar(data) && strcmp(data,'debug')
		unittests()
		return
	end

	Pfeas = [];
	Pinfeas = [];
	sizePopFeasible = configPrm.sizePopFeasible;
	sizePopInfeasible = configPrm.sizePopInfeasible;

	[nChunklets,chunklets] = generate_chunklets(conGraph);

	opts = statset('MaxIter',configPrm.maxKMSIter);

	nClusters = generate_nclusters(sizePopFeasible,configPrm);
	nTry = 0;
	while length(Pfeas) < sizePopFeasible && configPrm.maxInitTries > nTry
		nTry = nTry + 1;
		k = nClusters(length(Pfeas)+1);
		classOfClusters = spread_chunklets_among_clusters(nChunklets, k);
		seeds = sample_initial_seeds_off_chunklets(chunklets, classOfClusters);
		initMatrix = data(seeds,:);
		individual = initialize_gmm(data, k, initMatrix, opts);
		individual.classOfCluster = classOfClusters;
		[Pfeas Pinfeas] = insert_individual_correct_pool(individual, conGraph, Pfeas, Pinfeas,data);
	end

	nTry = 0;
	while length(Pinfeas) < sizePopInfeasible && configPrm.maxInitTries > nTry
		nTry = nTry + 1;
		k = nClusters(length(Pinfeas)+1);
		initMatrix = data(randsample(size(data,1),k),:);
		individual = initialize_gmm(data, k, initMatrix, opts);
		individual.classOfCluster = randi([1 configPrm.minClusters], k) ;

		[Pfeas Pinfeas] = insert_individual_correct_pool(individual, conGraph, Pfeas, Pinfeas,data);
	end

	Pfeas = Pfeas(1:min(length(Pfeas),sizePopFeasible));
	Pinfeas = Pinfeas(1:min(length(Pinfeas),sizePopInfeasible));
end


function [nChunklets outChunklets] = generate_chunklets(conGraph)
	MLs = conGraph>0;
	[nChunklets chunklets] = graphconncomp(MLs, 'Directed', false);
	validLabel = 1;
	outChunklets = chunklets;
	for i=1:nChunklets
		idx = find(chunklets == i);
		if numel(idx) == 1
			nChunklets = nChunklets - 1;
			outChunklets(idx) = 0;
		else
			outChunklets(idx) = validLabel;
			validLabel = validLabel + 1;
		end
	end
end


function [Pfeas Pinfeas] = insert_individual_correct_pool(individual, conGraph, Pfeas, Pinfeas,data)
	if isFeasible(individual, conGraph, data)
		idx = length(Pfeas);
		if idx == 0
			Pfeas = individual;
		else
			Pfeas(idx) = individual;
		end
	else
		idx =  length(Pinfeas);
		if idx == 0
			Pinfeas = individual;
		else
			Pinfeas(idx) = individual;
		end
	end
end

function isF = isFeasible(individual, conGraph, data)
	[i j s] = find(conGraph);
	constraints = [ i j s ];
	idxVio = compute_penalty(individual, constraints, data);
	isF = idxVio == false;
end

function individual = initialize_gmm(data, k, init, opts)
	[idx, clusters] = kmeans(data, k, 'EmptyAction', 'singleton', 'Start', init, 'Options', opts);
	individual = gmmFromKmeans(idx, clusters, data);
	individual = one_em_step(individual, data);
end

function individual = one_em_step(individual, data)
	individual = refinement(individual, data, struct('maxEMIter',1));
end

function seeds = sample_initial_seeds_off_chunklets(chunklets, classOfClusters)
	seeds = NaN([1 length(classOfClusters)]);
	for k=1:length(classOfClusters)
		candidates = find(chunklets==classOfClusters(k));
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


function [classOfClusters] = spread_chunklets_among_clusters(nChunklets, k)
	clustersToMap = k;
	classOfClusters = [];
	while clustersToMap > 0
		nToSample = min([nChunklets, clustersToMap]);
		classOfClusters = [classOfClusters;randsample(nChunklets,nToSample)];
		clustersToMap = clustersToMap - nToSample;
	end
end

function individual = create_empty_individual
		individual = struct( 'mean', [], ...
					 'covariance', [], ...
					 'mixCoef', [], ...
					 'distance', [], ...
					 'determinant',NaN,  ...
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
	testGenerateChunklets
	testInitializeWithConstraints
	testSpreadChunklets
	testGenerateNClusters
	testGmmFromKmeans
	testSampleInitialSeedsOfChunklets
end

function testGenerateChunklets
	graph = sparse([1 6 4], [3 7 5], [1 1 -1], 7, 7);
	graph = graph + graph';
	[outNchunk chunklets] = generate_chunklets(graph);
	assertEqual(outNchunk, 2)
	assertEqual(chunklets, [1 0 1 0 0 2 2])
end
function testInitializeWithConstraints
	data = [ 1 2; 1 3; 2 4; 5 5; 9 10; 10 11; 12 12];
	conHard = sparse([1 3 4 5], [3 5 6 6], [1 -1 -1 1], 7, 7);
	conHard = conHard+conHard';
	configPrm = struct('minClusters',2,'maxClusters',3, 'maxKMSIter',3, ...
		                 'sizePopFeasible',2, 'sizePopInfeasible', 2,...
	                   'maxInitTries', 5);

	[outFeas outInfeas] = initialize_with_constraints(data, conHard, configPrm)
	assertTrue(all(([outFeas(:).nClusters] >=2) & ([outFeas(:).nClusters] <= 3)),...
		'Numero de clusters invalido')
	for i=1:length(outFeas)
		[~,k1] = min(outFeas(i).distance(1,:));
		[~,k3] = min(outFeas(i).distance(3,:));
		[~,k4] = min(outFeas(i).distance(4,:));
		[~,k6] = min(outFeas(i).distance(6,:));
		assertEqual(outFeas(i).classOfCluster(k1), outFeas(i).classOfCluster(k3),...
			'Restrição ML violada em solucão factível')
		assertTrue(outFeas(i).classOfCluster(k4) ~= outFeas(i).classOfCluster(k6),...
			'Restrição CL violada em solucão factível')
	end

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
