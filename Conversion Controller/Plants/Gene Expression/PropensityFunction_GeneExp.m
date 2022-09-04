function Prop = PropensityFunction_GeneExp(x, Parameters)
% Propensity Function for Gene Expression
% 	 Species: 		 X = [X_1; X_2] 
% 	 Reactions: 	R1:         X_1            	--> 	0                   [gamma_1*X_1]
%                   R2:         X_2            	--> 	0                   [gamma_2*X_2]
%                   R3:         X_1            	--> 	X_1 + X_2           [k_1*X_1]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        gamma_2*X_2; ...
        k_1 * X_1; ...
       ];

end