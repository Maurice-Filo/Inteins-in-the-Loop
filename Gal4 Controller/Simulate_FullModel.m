%% Clear Workspace
close all
clear
clc
Save_Flag = 0;

%% Select Plant: {'Birth-Death', 'Gene Expression', 'Six Node'}
PlantChoice = 'Gene Expression';
switch PlantChoice
    % Birth-Death
    case 'Birth-Death'
        S_Plant = StoichiometryMatrix_BD();
        Prop_Plant = @PropensityFunction_BD;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.L = 1;
        Figure1_Name = 'FullModel_BD.pdf';
        
   	% Gene Expression
    case 'Gene Expression'
        S_Plant = StoichiometryMatrix_GeneExp();
        Prop_Plant = @PropensityFunction_GeneExp;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.gamma_2 = 1;
        Parameters_Plant.k_1 = 1;
        Parameters_Plant.L = 2;
        Figure1_Name = 'FullModel_GeneExp.pdf';
        
  	% Six Node Network
    case 'Six Node'
        S_Plant = StoichiometryMatrix_SixNode();
        Prop_Plant = @PropensityFunction_SixNode;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.gamma_2 = 0.1;
        Parameters_Plant.gamma_3 = 0.1;
        Parameters_Plant.gamma_4 = 0.1;
        Parameters_Plant.gamma_5 = 0.1;
        Parameters_Plant.gamma_6 = 0.1;
        Parameters_Plant.k_1 = 0.1;
        Parameters_Plant.k_2 = 0.1;
        Parameters_Plant.k_3 = 0.1;
        Parameters_Plant.k_4 = 0.1;
        Parameters_Plant.k_5 = 0.1;
        Parameters_Plant.gamma_F = 0.3;
        Parameters_Plant.kappa_F = 1;
        Parameters_Plant.L = 6;
        Figure1_Name = 'FullModel_SixNode.pdf';
end

%% Select Controller
[Parameters_Controller, S_Controller, Prop_Controller] = Design_Gal4Controller();
    epsilon = 0.01;
    Parameters_Controller.mu = [10; 0; 0; 0; 0; 0];
    Parameters_Controller.theta = [0; 1; 0; 0; 0; 0];
    Parameters_Controller.eta = 10;
    Parameters_Controller.a = [1; 1; 1] / epsilon;
    Parameters_Controller.d = [1; 1; 1] / epsilon;
    Parameters_Controller.c_1 = [] / epsilon;
    Parameters_Controller.c_2 = [] / epsilon;
    Parameters_Controller.delta_0 = 1;
    Parameters_Controller.delta = 0;
    Parameters_Controller.Parameters_Actuation.k = 0.01;

%% Simulation Settings
tf = 500;
N_t = 1000;
Solver = 'ODE23s';
DisturbanceFactor = 2;
t_Disturbance = tf/2;
Disturbed_Parameter = 'gamma_1';

%% Initial Conditions
IC_Full = zeros(Parameters_Controller.M + Parameters_Plant.L, 1);

%% Disturbance 
Parameters_DisturbedPlant = Parameters_Plant;
Parameters_DisturbedPlant.(Disturbed_Parameter) = Parameters_Plant.(Disturbed_Parameter) * DisturbanceFactor;

%% Closed-Loop Network for Full Controller Model
S_Mutual = zeros(Parameters_Plant.L, size(S_Controller,2));
S_Mutual(1,1:2) = [1, -1];
S_ClosedLoop = [S_Plant,                                            S_Mutual; ...
                zeros(Parameters_Controller.M, size(S_Plant,2)),  	S_Controller];
Parameters_ClosedLoop.Plant = Parameters_Plant;
Parameters_ClosedLoop.Controller = Parameters_Controller;
Prop_ClosedLoop = @(X, Parameters_ClosedLoop) ...
                    ([ Prop_Plant(X(1:Parameters_ClosedLoop.Plant.L), Parameters_ClosedLoop.Plant); ...
                       Prop_Controller(X(1), X(Parameters_ClosedLoop.Plant.L), X(Parameters_ClosedLoop.Plant.L+1:Parameters_ClosedLoop.Plant.L+Parameters_ClosedLoop.Controller.M), Parameters_ClosedLoop.Controller); ...
                    ]);
Parameters_DisturbedClosedLoop.Plant = Parameters_DisturbedPlant;
Parameters_DisturbedClosedLoop.Controller = Parameters_Controller;

%% Check RPA Conditions
q = Parameters_Controller.q_tot(:,1) - Parameters_Controller.q_tot(:,2);
if  all(q' * [Parameters_Controller.S_Q, Parameters_Controller.S_C, -Parameters_Controller.S_C, Parameters_Controller.S_B, -Parameters_Controller.S_B, Parameters_Controller.S_D] == 0)
    disp(['RPA conditions satisfied! (Setpoint = ', char(sym((-q'*Parameters_Controller.mu) / (q'*Parameters_Controller.theta))), ').']);
else
    disp('RPA conditions violated!');
end

%% Check Model Reduction Conditions
q_P = Parameters_Controller.q_tot(:,1);
q_N = Parameters_Controller.q_tot(:,2);
q_0 = Parameters_Controller.q_tot(:,3);
q_PN = [q_P, q_N];
q_tot = Parameters_Controller.q_tot;
m = size(Parameters_Controller.Q_1,1);
p = size(Parameters_Controller.B_1,1);
Condition11 = q_tot' * [Parameters_Controller.S_B, Parameters_Controller.S_C];
Condition12 = q_PN'* Parameters_Controller.S_D;
Condition21 = sign(q_P*q_N') - Parameters_Controller.Q_1'*Parameters_Controller.Q_2;
Condition22 = sign(q_P) - sign(Parameters_Controller.Q_1'*ones(m,1));
Condition23 = sign(q_N) - sign(Parameters_Controller.Q_2'*ones(m,1));
Condition31 = min(Parameters_Controller.Q_1*q_P, Parameters_Controller.Q_2*q_N) + Parameters_Controller.S_Q'*q_P;
Condition32 = min(Parameters_Controller.Q_1*q_P, Parameters_Controller.Q_2*q_N) + Parameters_Controller.S_Q'*q_N;
Condition33 = min(Parameters_Controller.Q_1*q_P, Parameters_Controller.Q_2*q_N) - Parameters_Controller.S_Q'*q_0;
Condition4 = max(Parameters_Controller.Q_1*q_P, Parameters_Controller.Q_2*q_N) - Parameters_Controller.v;
Condition5 = rank(Parameters_Controller.S_B) - size(Parameters_Controller.S_B,2);
Condition6 = rank([Parameters_Controller.S_B, Parameters_Controller.S_C]) - rank(Parameters_Controller.S_B) - rank(Parameters_Controller.S_C);
Condition7 = p + rank(Parameters_Controller.S_C) - (Parameters_Controller.M-3);
if ~any(Condition11, 'all') && ~any(Condition12, 'all') && ...
   ~any(Condition21, 'all') && ~any(Condition22, 'all') && ~any(Condition23, 'all') && ...
   ~any(Condition31, 'all') && ~any(Condition32, 'all') && ~any(Condition33, 'all') && ...
   ~any(Condition4, 'all') && ~any(Condition5, 'all') && ~any(Condition6, 'all') && ~any(Condition7, 'all')
    disp('Model reduction conditions satisfied!');
end

%% Simulating the Full Model
[t_Full_1, x_Full_1] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_ClosedLoop, IC_Full, t_Disturbance, N_t, Solver);
[t_Full_2, x_Full_2] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_DisturbedClosedLoop, x_Full_1(:,end), tf - t_Disturbance, N_t, Solver);
t_Full = [t_Full_1, t_Disturbance + t_Full_2(2:end)];
x_Full = [x_Full_1, x_Full_2(:,2:end)];

%% Figure Settings
Scale = 4;
Figure_Width = 10 * Scale;
Figure_Height = 10 * Scale;
FontSize = 5 * Scale;
FontSize_Small = 3 * Scale;
FontSize_Large = 6 * Scale;
LineWidth = 0.65 * Scale;
LineWidth_Thick = 0.8 * Scale;
LineWidth_Thin = 0.1 * Scale;
MarkerSize = 5 * Scale;
Opacity = 0.5;
Colors = lines(20);
    
%% Set Figure 1
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Simulation Plots
Species = cell(Parameters_Plant.L + Parameters_Controller.M,1);
for i = 1 : Parameters_Plant.L
    Species{i} = ['X_', num2str(i)];
end
for i = 1 : Parameters_Controller.M
    Species{i+Parameters_Plant.L} = ['Z_', num2str(i)];
end
Handle_Tile = tiledlayout(Handle_Figure1, ceil((Parameters_Plant.L + Parameters_Controller.M)/2), 2);
    Handle_Tile.TileSpacing = 'compact';
    Handle_Tile.Padding = 'compact';
    title(Handle_Tile, 'Response', 'FontSize', FontSize, 'FontWeight', 'bold');
    xlabel(Handle_Tile, 'Time', 'FontSize', FontSize, 'FontWeight', 'bold');
    ylabel(Handle_Tile, 'Concentration', 'FontSize', FontSize, 'FontWeight', 'bold');
for i = 1 : Parameters_Plant.L + Parameters_Controller.M
    Handle_Axis = nexttile;
    plot(t_Full, x_Full(i, :), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick);
    title(['$\bf{', Species{i}, '}$']);
    Handle_Axis.XLim = [0, tf];
    Handle_Axis.Title.Interpreter = 'latex';
    grid(Handle_Axis, 'on');
    Handle_Axis.FontSize = FontSize;
end

%% Save Figure
if Save_Flag == 1
    exportgraphics(Handle_Tile, Figure1_Name, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
