function conGraph = generate_constraint_graph(constraints, nObjects)
%GENERATE_CONSTRAINT_GRAPH Given a constraint matrix Cx3 transform it in its graph representations
	conGraph= sparse(constraints(:,1), constraints(:,2), constraints(:,3), nObjects, nObjects );
	%make it symmetric
	conGraph = conGraph + conGraph';
