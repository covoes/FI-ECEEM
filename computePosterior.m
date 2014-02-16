function [posterior gauss] = computePosterior( individual, data )
%COMPUTEPOSTERIOR Computer the posterior probability for each individual to each cluster
%
%In the second argument the gaussian pdfs for the objects are returned in
%a nObjects x nClusters matrix

[nObjects nFeatures] = size(data);
nClusters = individual.nClusters;
mixCoef = individual.mixCoef(1:nClusters);

posterior = NaN([nObjects nClusters]);

gauss = NaN([nObjects nClusters]);

const = (2*pi)^(nFeatures/2);
for k=1:nClusters
	gauss(:,k) = (1/(const* (individual.determinant(k))^0.5)) .* exp(-0.5.*individual.distance(:,k));
	posterior(:,k) = mixCoef(k) .* gauss(:,k);
end

%normalizing posteriors
sums = sum(posterior, 2);
posterior = bsxfun(@rdivide, posterior, sums);

end
