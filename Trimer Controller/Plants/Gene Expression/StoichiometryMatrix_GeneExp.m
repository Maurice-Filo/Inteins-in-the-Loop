function S = StoichiometryMatrix_GeneExp()
% Stoichiometry Matrix for Gene Expression
% 	 Species: 		 X = [X_1; X_2] 
% 	 Reactions: 	R1:         X_1            	--> 	0                   [gamma_1*X_1]
%                   R2:         X_2            	--> 	0                   [gamma_2*X_2]
%                   R3:         X_1            	--> 	X_1 + X_2           [k_1*X_1]

S = [ ...
    -1     0     0; ...          
     0    -1     1; ...          
    ];
end
