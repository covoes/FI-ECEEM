function [idxVio totPenalty] = compute_penalty(P, constraints, data)
%COMPUTE_PENALTY Given a population, compute a penalty score for each individual and return in 
%    idxVio a boolean vector of individuals with constraint violations

idxVio = logical(zeros([1 length(p)]))
totPenalty = zeros([1 length(p)])

for i=1:length(P)
	for c=1:size(constraints,1)
		idx1 = constraints(c,1)
		idx2 = constraints(c,2)
		[~,gauss] = computePosterior(P(i), data)
		%pk1 holds the responsability of gaussian k1 generate idx1
		%pk2 holds the analogue for idx2 [for k2]
		[pk1,k1] = max(gauss(idx1,:))
		[pk2,k2] = max(gauss(idx2,:))
		classK1 = P(i).classOf(k1)
		classK2 = P(i).classOf(k2)
		penalty = 0
		if constraints(c,3) == 1 && classK1 != classK2
			%ML constraint being violated
			penalty = (1-gauss(idx2, k1)) + (1-gauss(idx1, k2))
		else if constraints(c,3) == -1 && classK1 == classK2
			%CL constraint being violated
      penalty = gauss(idx1, k1) +  gauss(idx2, k1)
		end
    totPenalty(i) = totPenalty(i) + penalty
	end
end
idxVio(totPenalty>0) = true
