function [u_P, u_N, z_tilde] = ReducedControlAction_intraDDController(~, z_p, ~, ~, Parameters)
%% Extract Parameters
kappa = Parameters.kappa;
k = Parameters.Parameters_Actuation.k;

%% Compute Original State Variables
z_tilde = (1/8) * ( 4*z_p + kappa - sqrt(8*kappa*z_p + kappa^2) );
    
%% Compute Control Action u
u_P = k*z_tilde;  
u_N = 0;
end