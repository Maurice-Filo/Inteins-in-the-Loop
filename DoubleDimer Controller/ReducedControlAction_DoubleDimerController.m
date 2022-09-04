function [u_P, u_N, z_tilde] = ReducedControlAction_DoubleDimerController(~, z_p, z_n, z_0, Parameters)
%% Extract Parameters
kappa_1 = Parameters.kappa(1);
kappa_2 = Parameters.kappa(2);
kappa_3 = Parameters.kappa(3);
k = Parameters.Parameters_Actuation.k;
kappa_u = Parameters.Parameters_Actuation.kappa_u;
kappa_u_prime = Parameters.Parameters_Actuation.kappa_u_prime;
    
%% Compute Original State Variables
z_tilde_1 = (1/8) * (4*z_p + kappa_1 - sqrt(8*kappa_1*z_p + kappa_1^2) );
z_tilde_2 = (1/8) * (4*z_n + kappa_2 - sqrt(8*kappa_2*z_n + kappa_2^2) );
z_tilde_3 = (1/8) * (4*z_0 + kappa_3 - sqrt(8*kappa_3*z_0 + kappa_3^2) );
z_tilde = [z_tilde_1; z_tilde_2; z_tilde_3];
    
%% Compute Control Action
u_P = (k*z_tilde_1) / (1 + z_tilde_2/kappa_u + z_tilde_3/kappa_u_prime);
u_N = 0;
end

