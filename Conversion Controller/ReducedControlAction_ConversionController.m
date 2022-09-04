function [u_P, u_N, z_tilde] = ReducedControlAction_ConversionController(~, z_p, z_n, ~, Parameters)
%% Extract Parameters
rho_1 = Parameters.rho(1);
rho_2 = Parameters.rho(2);
k = Parameters.Parameters_Actuation.k;
gamma = Parameters.Parameters_Actuation.gamma;
    
%% Compute Original State Variables
z_tilde_1 = z_p / (1 + rho_1);
z_tilde_2 = z_n / (1 + rho_2);
z_tilde = [z_tilde_1; z_tilde_2];
    
%% Compute Control Action u
u_P = k*z_tilde_1;  
u_N = gamma * z_tilde_2;
end

