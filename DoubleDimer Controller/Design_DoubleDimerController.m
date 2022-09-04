function [Parameters_Controller, S, Prop] = Design_DoubleDimerController()
%% Number of Controller Species and Charges
M = 6; 
q_P  = [1, 0, 0, 2, 0, 0]';
q_N  = [0, 1, 0, 0, 2, 0]';
q_0   = [0, 0, 1, 0, 0, 2]';
q_tot = [q_P, q_N, q_0];

%% Production Reactions;
mu = [10; 0; 0; 0; 0; 0];
theta = [0; 1; 0; 0; 0; 0];

%% Irreversible Sequestration Reactions: Z_i + Z_j -> Z_k + ...
SequestrationRxns = { ...
                        'Z_1 + Z_2 -> Z_3'; ...
                        'Z_4 + Z_5 -> Z_6'; ...
                        'Z_4 + Z_2 -> Z_1 + Z_3'; ...
                        'Z_1 + Z_5 -> Z_2 + Z_3'; ...
                    };  
v = [1; 2; 2; 2];
eta = 10; 
[S_Q, Q_1, Q_2, Q_3] = Sequestration(SequestrationRxns, M);

%% Reversible Binding Reactions: Z_i + Z_j <-> Z_k
BindingRxns = { ...
                'Z_1 + Z_1 <-> Z_4'; ...
                'Z_2 + Z_2 <-> Z_5'; ...
                'Z_3 + Z_3 <-> Z_6'; ...
              };
a = [];
d = [];
[S_B, B_1, B_2, B_3] = Binding(BindingRxns, M);

%% Reversible Conversion Reactions: Z_i <-> Z_j
ConversionRxns =    { ...
                    }; 
c_1 = [];
c_2 = [];
[S_C, C_1, C_2] = Conversion(ConversionRxns, M);

%% Degradation Reactions: Z_i -> 0
DegradationRxns =   { ...
                        'Z_3 -> 0'; ...
                        'Z_6 -> 0'; ...
                    };
delta_0 = 1;
[S_D, Q_0] = Degradation(DegradationRxns, M);

%% Actuation Functions
Parameters_Actuation.k = 0.01; 
Parameters_Actuation.kappa_u = 1;
Parameters_Actuation.kappa_u_prime = 1;
h_P = @(x_L, z, ActuationParameters) (ActuationParameters.k * z(4) / (1 + z(5)/Parameters_Actuation.kappa_u + z(6)/Parameters_Actuation.kappa_u_prime));
h_N = @(x_L, z, ActuationParameters) (0);

%% Dilution Rate
delta = 0;

%% Store Controller Parameters
Parameters_Controller.M = M;
Parameters_Controller.q_tot = q_tot;
Parameters_Controller.mu = mu;
Parameters_Controller.theta = theta;
Parameters_Controller.v = v;
Parameters_Controller.eta = eta;
Parameters_Controller.a = a;
Parameters_Controller.d = d;
Parameters_Controller.c_1 = c_1;
Parameters_Controller.c_2 = c_2;
Parameters_Controller.delta_0 = delta_0;
Parameters_Controller.delta = delta;
Parameters_Controller.Parameters_Actuation = Parameters_Actuation;
Parameters_Controller.h_P = h_P;
Parameters_Controller.h_N = h_N;
Parameters_Controller.Q_1 = Q_1;
Parameters_Controller.Q_2 = Q_2;
Parameters_Controller.Q_3 = Q_3;
Parameters_Controller.S_Q = S_Q;
Parameters_Controller.B_1 = B_1;
Parameters_Controller.B_2 = B_2;
Parameters_Controller.B_3 = B_3;
Parameters_Controller.S_B = S_B;
Parameters_Controller.C_1 = C_1;
Parameters_Controller.C_2 = C_2;
Parameters_Controller.S_C = S_C;
Parameters_Controller.Q_0 = Q_0;
Parameters_Controller.S_D = S_D;

%% Controller Stoichiometry Matrix
S = [zeros(M,2), eye(M), S_Q, S_C, -S_C, S_B, -S_B, S_D, -eye(M)];

%% Controller Propensity Function
Prop = @(x_1, x_L, z, Parameters_Controller) ...
        ([  Parameters_Controller.h_P(x_L, z, Parameters_Controller.Parameters_Actuation); ...                                           	% Positive Actuation
            Parameters_Controller.h_N(x_L, z, Parameters_Controller.Parameters_Actuation)*x_1; ...                                        	% Negative Actuation
         	Parameters_Controller.mu + Parameters_Controller.theta*x_L; ...                                                               	% Setpoint/Sensing
            (Parameters_Controller.eta*Parameters_Controller.v) .* (Parameters_Controller.Q_1*z) .* (Parameters_Controller.Q_2*z); ...    	% Sequestration
          	Parameters_Controller.c_1 .* (Parameters_Controller.C_1*z); ...                                                                 % Forward Conversion         
         	Parameters_Controller.c_2 .* (Parameters_Controller.C_2*z); ...                                                                 % Backward Conversion  
          	Parameters_Controller.a .* (Parameters_Controller.B_1*z) .* (Parameters_Controller.B_2*z); ...                                  % Forward Binding
         	Parameters_Controller.d .* (Parameters_Controller.B_3*z); ...                                                                   % Backward Binding              
          	Parameters_Controller.delta_0 * (Parameters_Controller.Q_0*z); ...                                                              % Degradation
          	Parameters_Controller.delta*z; ...                                                                                              % Dilution
       	]);
    
end

