function Prop = PropensityFunction_SixNode(x, Parameters)
% Propensity Function for the Six Node Network
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

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
gamma_4 = Parameters.gamma_4;
gamma_5 = Parameters.gamma_5;
gamma_6 = Parameters.gamma_6;
k_1 = Parameters.k_1;
k_2 = Parameters.k_2;
k_3 = Parameters.k_3;
k_4 = Parameters.k_4;
k_5 = Parameters.k_5;
gamma_F = Parameters.gamma_F;
kappa_F = Parameters.kappa_F;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);
X_4 = x(4);
X_5 = x(5);
X_6 = x(6);

%% Propensities
Prop = [ ...    
        gamma_1*X_1; ...
        gamma_2*X_2; ...
        gamma_3*X_3; ...
        gamma_4*X_4; ...
        gamma_5*X_5; ...
        gamma_6*X_6; ...
        k_1*X_1; ...
        k_2*X_2; ...
        k_3*X_3; ...
        k_4*X_4; ...
        k_5*X_5; ...
        gamma_F*X_6*X_2/(X_2 + kappa_F); ...
       ];

end