function [u_P, u_N, z_tilde] = ReducedControlAction_Gal4Controller(~, z_p, ~, z_0, Parameters)
%% Extract Parameters
kappa_1 = Parameters.kappa(1);
kappa_2 = Parameters.kappa(2);
kappa_3 = Parameters.kappa(3);
k = Parameters.Parameters_Actuation.k;

%% Polynomial Coefficients: c_4*y_1^4 + c_3*y_1^3 + c_2*y_1^2 + c_1*y_1 + c_0
c_4 = 8*kappa_2^2 - 2*kappa_1*kappa_3;
c_3 = 8*kappa_2^2*sqrt(kappa_3) - kappa_1*kappa_3^(3/2) - 2*kappa_1*kappa_2*sqrt(kappa_3);
c_2 = 2*kappa_2^2*kappa_3 - 8*kappa_2^2*z_0 - kappa_1*kappa_2*kappa_3 + kappa_1*kappa_3*z_0 - kappa_1*kappa_3*z_p;
c_1 = kappa_1*kappa_2*sqrt(kappa_3)*z_0 - 4*kappa_2^2*sqrt(kappa_3)*z_0;
c_0 = 2*kappa_2^2*z_0^2;
Coefficients = [c_4, c_3, c_2, c_1, c_0];
    
%% Solve Polynomial in Transformed Variables [y_1, y_2]
if z_p < eps
    y_2 = 0; y_1 = (1/4) * (-sqrt(kappa_3) + sqrt(kappa_3 + 8*z_0) );
elseif z_0 < eps
    y_1 = 0; y_2 = (1/4) * (-sqrt(kappa_1) + sqrt(kappa_1 + 8*z_p) );
else
    y_1 = roots(Coefficients);
    y_1 = y_1(real(y_1)>=0 & imag(y_1)==0);
    y_2 = (kappa_2 / sqrt(kappa_1*kappa_3) ) * ( (z_0 - 2*y_1.^2 - sqrt(kappa_3)*y_1) ./ y_1 );
    y_1 = y_1(real(y_2)>=0 & imag(y_2)==0);
    y_2 = y_2(real(y_2)>=0 & imag(y_2)==0);
end

if length(y_1) ~= 1 || length(y_2) ~= 1
    keyboard
end
    
%% Compute Original State Variables
z_tilde_1 = y_1^2;
z_tilde_2 = y_2^2;
z_tilde_3 = sqrt(z_tilde_1*z_tilde_2);
z_tilde = [z_tilde_1; z_tilde_2; z_tilde_3];
    
%% Compute Control Action u
u_P = k*z_tilde_2;  
u_N = 0;
end