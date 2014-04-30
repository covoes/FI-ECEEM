function [fitness] = fitnessFunc( FName, NlogL, numObjects, numClusters, numFeatures )
%Compute the fitness of a gmdistribution object.
%
%Possible values for FName are: 'mdl', 'fig','fig2'
%If adding a new function make sure that its in a minimization way
%i.e. the intention is to minimize the function

switch( lower(FName) )

	case 'mdl'
		penalMDL = (numClusters ... %mixing coefficients
			 + numClusters * numFeatures ... %means
		     + numClusters * numFeatures * (numFeatures+1)/2)/2 ... %covariance matrix
			 * log(numObjects); %data
		fitness = NlogL + penalMDL;
		%sprintf('%.5f == %.5f\n', fitness*2, GMDistObj.BIC)
		%fitness = GMDistObj.BIC;

	otherwise
		error( 'Function %s not  implemented', FName );
	end
end
