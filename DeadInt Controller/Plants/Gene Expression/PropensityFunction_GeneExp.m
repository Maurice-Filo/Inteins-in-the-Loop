function Propensity = PropensityFunction_GeneExp(x, Parameters)
% Propensity Function for Gene Expression Network
% 	 Species:           X = [X_1; X_2]
% 	 Reactions:         R1:		 X_1					-->         X_1 +  X_2				[k_1*X_1]
%                       R2:		 X_1					-->         phi                     [gamma_p_1*X_1]
%                       R3:		 X_2					-->         phi                     [gamma_p_2*X_2]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_p_1 = Parameters.gamma_p_1;
gamma_p_2 = Parameters.gamma_p_2;

%% Extract Variables
X_1 = x(1);
X_2 = x(2);

%% Construct Propensity Function
Propensity = [  k_1 * X_1; ...
                gamma_p_1 * X_1; ...
                gamma_p_2 * X_2 ...
              ];
end

