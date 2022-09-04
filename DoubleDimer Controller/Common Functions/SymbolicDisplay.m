function [RHS] = SymbolicDisplay(L, M, Parameters_CL, S_CL, Prop_CL)
%% Define Symbolic Variables
x = (sym('x_', [1, L], 'positive'))';
z = (sym('z_', [1, M], 'positive'))';
Variables = [x; z];

%% Define Symbolic Parameters
% Plant Parameters
PlantParams = fieldnames(Parameters_CL.Plant);
for i = 1 : length(PlantParams)
    if ~strcmp(PlantParams{i}, 'L')
        eval([PlantParams{i}, ' = sym(''', PlantParams{i}, ''', ''positive'');']); 
        eval(['Parameters.Plant.', PlantParams{i}, ' = ', PlantParams{i}, ';']);
    end
end
% Actuation Parameters
ActParams = fieldnames(Parameters_CL.Controller.Parameters_Actuation);
for i = 1 : length(ActParams)
    eval([ActParams{i}, ' = sym(''', ActParams{i}, ''', ''positive'');']); 
    eval(['Parameters.Controller.Parameters_Actuation.', ActParams{i}, ' = ', ActParams{i}, ';']);
end
Parameters.Controller.h_P = Parameters_CL.Controller.h_P;
Parameters.Controller.h_N = Parameters_CL.Controller.h_N;
% Binding Parameters
p = length(Parameters_CL.Controller.a);
Parameters.Controller.a = (sym('a', [1 p], 'positive'))';
Parameters.Controller.d = (sym('d', [1 p], 'positive'))';
% Sequestration Parameters
Parameters.Controller.eta = (sym('eta', 'positive'))';
Parameters.Controller.v = Parameters_CL.Controller.v;
% Conversion Parameters
p_0 = length(Parameters_CL.Controller.c_1);
Parameters.Controller.c_1 = (sym('c_1_', [1 p_0], 'positive'))';
Parameters.Controller.c_2 = (sym('c_2_', [1 p_0], 'positive'))';
% Conversion Parameters
Parameters.Controller.delta_0 = (sym('delta_0', 'positive'))';
% Dilution Parameter
syms delta positive;
Parameters.Controller.delta = delta;
% Production Parameters
mu = (sym('mu_', [1, M], 'positive'))';
theta = (sym('theta_', [1, M], 'positive'))';
for i = 1 : M
    if Parameters_CL.Controller.mu(i) == 0 
        mu(i) = 0;
    end
    if Parameters_CL.Controller.theta(i) == 0 
        theta(i) = 0;
    end
end
Parameters.Controller.mu = mu;
Parameters.Controller.theta = theta;
% Matrices
Parameters.Controller.Q_1 = Parameters_CL.Controller.Q_1;
Parameters.Controller.Q_2 = Parameters_CL.Controller.Q_2;
Parameters.Controller.C_1 = Parameters_CL.Controller.C_1;
Parameters.Controller.C_2 = Parameters_CL.Controller.C_2;
Parameters.Controller.B_1 = Parameters_CL.Controller.B_1;
Parameters.Controller.B_2 = Parameters_CL.Controller.B_2;
Parameters.Controller.B_3 = Parameters_CL.Controller.B_3;
Parameters.Controller.Q_0 = Parameters_CL.Controller.Q_0;

%% Dimensions
Parameters.Controller.M = M;
Parameters.Plant.L = L;

%% ODEs
RHS = S_CL * Prop_CL(Variables, Parameters);
RHS_Latex = cell(length(RHS),1);
for i = 1 : length(RHS)
	RHS_Latex{i} = ['$ \dot ', char(Variables(i)), ' = ', latex(RHS(i)), '$'];

end
figure();
h = annotation(gcf,'textbox', [0, 0, 1, 1], 'String', RHS_Latex, 'interpreter', 'latex', 'FontSize', 30);
end


