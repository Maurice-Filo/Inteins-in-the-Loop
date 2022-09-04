function S = StoichiometryMatrix_SixNode()
% Stoichiometry Matrix for the Six Node Network
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6]
% 	 Reactions: 	R1:         X_1				--> 	0                   [gamma_1*X_1]
% 				    R2:         X_2				--> 	0                   [gamma_2*X_2]
% 				    R3:         X_3				--> 	0                   [gamma_3*X_3]
% 				    R4:         X_4				--> 	0                   [gamma_4*X_4]
% 				    R5:         X_5				--> 	0                   [gamma_5*X_5]
% 				    R6:         X_6				--> 	0                   [gamma_6*X_6]
% 				    R7:         X_1				--> 	X_1 +  X_2          [k_1*X_1]
% 				    R8:         X_2				--> 	X_2 +  X_3          [k_2*X_2]
% 				    R9:         X_3				--> 	X_3 +  X_4          [k_3*X_3]
% 				    R10:		X_4				--> 	X_4 +  X_5          [k_4*X_4]
% 				    R11:		X_5				--> 	X_5 +  X_6          [k_5*X_5]
% 				    R12:		X_2				--> 	0                   [gamma_F*X_6*X_2/(X_2 + kappa_F)]

S = [ ...
    -1     0	 0	   0	 0	   0	 0	   0	 0	   0	 0	   0; ...
     0    -1	 0	   0	 0	   0	 1	   0	 0     0	 0	  -1; ...
     0     0	-1	   0	 0	   0	 0	   1	 0	   0	 0	   0; ...
     0	   0	 0	  -1	 0	   0	 0	   0	 1	   0	 0	   0; ...
     0	   0	 0	   0	-1	   0	 0	   0	 0	   1	 0	   0; ...
     0	   0	 0	   0	 0	  -1	 0	   0	 0	   0	 1	   0; ...
    ];


end