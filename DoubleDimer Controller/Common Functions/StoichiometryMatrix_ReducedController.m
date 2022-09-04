function S = StoichiometryMatrix_ReducedController()
% Stoichiometry Matrix for the Reduced Controller
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


S = [ ...          
     0     0     1     0     0    -1    -1     0     0; ...
     0     0     0     1     0    -1     0    -1     0; ...
     0     0     0     0     1     1     0     0    -1; ...
    ];
end
