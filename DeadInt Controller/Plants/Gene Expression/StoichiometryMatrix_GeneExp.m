function S = StoichiometryMatrix_GeneExp()
% Stoichiometry Matrix for Gene Expression Network
% 	 Species:           X = [X_1; X_2]
% 	 Reactions:         R1:		 X_1					-->         X_1 +  X_2				[k_1*X_1]
%                       R2:		 X_1					-->         phi                     [gamma_p_1*X_1]
%                       R3:		 X_2					-->         phi                     [gamma_p_2*X_2]
S = [	0,	   -1,		0; ...
		1,		0,	   -1; ...
	];

end

