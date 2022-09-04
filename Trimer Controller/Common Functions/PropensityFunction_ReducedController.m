function Prop = PropensityFunction_ReducedController(X_1, X_L, z, Parameters)
% Propensity Function for the Reduced Controller
% 	 Species: 		 X = [Z_P; Z_N; Z_0] 
% 	 Reactions: 	R1:         0               -->     X_1                 [u_P]
%                   R2:         X_1             -->     0                   [u_N*x_1]
%                   R3:         0               --> 	Z_P                 [mu_P + theta_P*x_L]
%                   R4:         0               -->     Z_N                 [mu_N + theta_N*X_L]
%                   R5:         0               -->     Z_0                 [mu_0 + theta_0*X_L + delta_0*q_PN0'*phi]
%                   R6:         Z_P + Z_N       -->     Z_0                 [eta*Z_P*Z_N]
%                   R7:         Z_P             -->     0                   [delta*Z_P]
%                   R8:         Z_N             -->     0                   [delta*Z_N]
%                   R9:         Z_0             -->     0                   [(delta + delta_0)*Z_0]

%% Extract Parameters
mu_P = Parameters.mu_P;
theta_P = Parameters.theta_P;
mu_N = Parameters.mu_N;
theta_N = Parameters.theta_N;
mu_0 = Parameters.mu_0;
theta_0 = Parameters.theta_0;
eta = Parameters.eta;
delta = Parameters.delta;
delta_0 = Parameters.delta_0;
ControlAction = Parameters.ControlAction;
q_tot = Parameters.q_tot;
W_1 = Parameters.W_1;
W_2 = Parameters.W_2;

%% Extract State Variables
Z_P = z(1);
Z_N = z(2);
Z_0 = z(3);

%% Propensities
[u_P, u_N, Z_tilde] = ControlAction(X_L, Z_P, Z_N, Z_0, Parameters);
Z_tot = [Z_P; Z_N; Z_0];
Term = delta_0 * ((q_tot(:,1) + q_tot(:,2)) .* q_tot(:,3))' * (W_1' * Z_tot + W_2' * Z_tilde);
Prop = [ ...
        u_P; ...  
        u_N*X_1; ...
        mu_P + theta_P*X_L; ...
        mu_N + theta_N*X_L; ...
        mu_0 + theta_0*X_L + Term; ...
       	eta*Z_P*Z_N; ...
        delta*Z_P; ...
        delta*Z_N; ...
     	(delta + delta_0)*Z_0; ...
       ];
end