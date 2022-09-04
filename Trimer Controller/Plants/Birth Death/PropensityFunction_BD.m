function Prop = PropensityFunction_BD(x, Parameters)
% Propensity Function for the Birth-Death Process
% 	 Species: 		 X = [X_1]
% 	 Reactions: 	R1:         X_1            	--> 	0                   [gamma_1*X_1]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;

%% Extract State Variables
X_1 = x(1);

%% Propensities
Prop = gamma_1*X_1;

end