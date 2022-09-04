%% Clear Workspace
% close all
clear
clc

%% Parameters
a = 10*[1; 1; 1; 1; 1; 1];
d = 10*[1; 1; 1; 1; 1; 1];
k = 1;
k_prime = 0.5;
kappa_u = 1;

%% Inputs
N = 30; M = 30;
z_p_vector = linspace(0, 60, N);
z_s_vector = linspace(0, 70, M);

%% Evaluate
IG = zeros(4,1);
z_5 = zeros(N,M);
z_6 = zeros(N,M);
z_7 = zeros(N,M);
z_8 = zeros(N,M);
for i = 1 : N
    z_p = z_p_vector(i);
    for j = 1 : M
        z_s = z_s_vector(j);
        [z_5(i,j), z_6(i,j), z_7(i,j), z_8(i,j)] = AlgebraicEqn(z_p, z_s, a, d, IG);
        IG = [z_5(i,j), z_6(i,j), z_7(i,j), z_8(i,j)]';
    end
end
phi_z = repmat(z_p_vector', 1, M) - (z_6 + 2*z_7 + z_8);
psi_z = 2*z_6 + 2*z_7 + z_8;
u_z = (k*z_7 + k_prime*z_6) ./ (1 + z_5 / kappa_u);

%% Plotting Results
figure
subplot(1,3,1)
surf(z_p_vector, z_s_vector, phi_z');
subplot(1,3,2)
surf(z_p_vector, z_s_vector, psi_z');
subplot(1,3,3)
surf(z_p_vector, z_s_vector, u_z');