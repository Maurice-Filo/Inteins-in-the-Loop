%% Clear Workspace
close all
clear
clc
Save_Flag = 1;

%% Select Plant: {'Birth-Death', 'Gene Expression', 'Six Node'}
PlantChoice = 'Six Node';
switch PlantChoice
    % Birth-Death
    case 'Birth-Death'
        S_Plant = StoichiometryMatrix_BD();
        Prop_Plant = @PropensityFunction_BD;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.L = 1;
        Figure1_Name = 'ModelReduction_Gal4_BD.pdf';
        % Swept Parameter
        ParameterName = 'gamma_1';
        ParameterValues = [0.1, 0.5, 1, 1.5];
        
   	% Gene Expression
    case 'Gene Expression'
        S_Plant = StoichiometryMatrix_GeneExp();
        Prop_Plant = @PropensityFunction_GeneExp;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.gamma_2 = 1;
        Parameters_Plant.k_1 = 1;
        Parameters_Plant.L = 2;
        Figure1_Name = 'ModelReduction_Gal4_GeneExp.pdf';
        % Swept Parameter
        ParameterName = 'gamma_2';
        ParameterValues = [0.1, 0.3, 0.5, 1];
        
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
        Figure1_Name = 'ModelReduction_Gal4_SixNode.pdf';
        % Swept Parameter
        ParameterName = 'gamma_6';
        ParameterValues = [0.1, 0.5, 1, 1.5];
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
    Parameters_Controller.delta = 0.005;
    Parameters_Controller.Parameters_Actuation.k = 0.05;
S_ReducedController = StoichiometryMatrix_ReducedController();
Prop_ReducedController = @PropensityFunction_ReducedController;

%% Select Reduced Control Action
Parameters_ReducedController.ControlAction = @ReducedControlAction_Gal4Controller;

%% Simulation Settings
tf = 500;
N_t = 1000;
Solver = 'ODE23s';
DisturbanceFactor = 2;
t_Disturbance = tf/2;
Disturbed_Parameter = 'gamma_1';

%% Initial Conditions
IC_Full = zeros(Parameters_Controller.M + Parameters_Plant.L, 1);
TransformationMatrix = [eye(Parameters_Plant.L), zeros(Parameters_Plant.L, Parameters_Controller.M); ...
                        zeros(3, Parameters_Plant.L), Parameters_Controller.q_tot'];
IC_Reduced = TransformationMatrix*IC_Full;

%% Disturbance 
Parameters_DisturbedPlant = Parameters_Plant;
Parameters_DisturbedPlant.(Disturbed_Parameter) = Parameters_Plant.(Disturbed_Parameter) * DisturbanceFactor;

%% Parameters of Reduced Model
Parameters_ReducedController.mu_P = Parameters_Controller.q_tot(:,1)' * Parameters_Controller.mu;
Parameters_ReducedController.theta_P = Parameters_Controller.q_tot(:,1)' * Parameters_Controller.theta;
Parameters_ReducedController.mu_N = Parameters_Controller.q_tot(:,2)' * Parameters_Controller.mu;
Parameters_ReducedController.theta_N = Parameters_Controller.q_tot(:,2)' * Parameters_Controller.theta;
Parameters_ReducedController.mu_0 = Parameters_Controller.q_tot(:,3)' * Parameters_Controller.mu;
Parameters_ReducedController.theta_0 = Parameters_Controller.q_tot(:,3)' * Parameters_Controller.theta;
Parameters_ReducedController.eta = Parameters_Controller.eta;
Parameters_ReducedController.delta_0 = Parameters_Controller.delta_0;
Parameters_ReducedController.delta = Parameters_Controller.delta;
C_bar = [eye(3), zeros(3,Parameters_Controller.M-3)];
C_tilde = [zeros(Parameters_Controller.M-3,3), eye(Parameters_Controller.M-3)];
Parameters_ReducedController.W_1 = (C_bar * Parameters_Controller.q_tot) \ C_bar;
Parameters_ReducedController.W_2 =  C_tilde * (eye(Parameters_Controller.M) -  Parameters_Controller.q_tot*Parameters_ReducedController.W_1);
Parameters_ReducedController.kappa = Parameters_Controller.d ./ Parameters_Controller.a;
Parameters_ReducedController.rho = Parameters_Controller.c_2 ./ Parameters_Controller.c_1;
Parameters_ReducedController.q_tot = Parameters_Controller.q_tot;
Parameters_ReducedController.Parameters_Actuation = Parameters_Controller.Parameters_Actuation;

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

%% Closed-Loop Network for Reduced Controller Model
S_Mutual_Reduced = zeros(Parameters_Plant.L, size(S_ReducedController,2));
S_Mutual_Reduced(1,1:2) = [1, -1];
S_ReducedClosedLoop = [S_Plant,                       	S_Mutual_Reduced; ...
                       zeros(3, size(S_Plant,2)),  	S_ReducedController];
Parameters_ReducedClosedLoop.Plant = Parameters_Plant;
Parameters_ReducedClosedLoop.Controller = Parameters_ReducedController;
Prop_ReducedClosedLoop = @(X, Parameters_ReducedClosedLoop) ...
                    ([ Prop_Plant(X(1:Parameters_ReducedClosedLoop.Plant.L), Parameters_ReducedClosedLoop.Plant); ...
                       Prop_ReducedController(X(1), X(Parameters_ReducedClosedLoop.Plant.L), X(Parameters_ReducedClosedLoop.Plant.L+1:Parameters_ReducedClosedLoop.Plant.L+3), Parameters_ReducedClosedLoop.Controller); ...
                    ]);
Parameters_DisturbedReducedClosedLoop.Plant = Parameters_DisturbedPlant;
Parameters_DisturbedReducedClosedLoop.Controller = Parameters_ReducedController;

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

%% Simulations
x_Full_Transformed = zeros(Parameters_Plant.L + 3, 2*N_t - 1, length(ParameterValues));
x_Reduced = zeros(Parameters_Plant.L + 3, 2*N_t - 1, length(ParameterValues));

for Counter = 1 : length(ParameterValues)
    Parameters_ClosedLoop.Plant.(ParameterName) = ParameterValues(Counter);
    Parameters_DisturbedClosedLoop.Plant.(ParameterName) = ParameterValues(Counter);
    Parameters_ReducedClosedLoop.Plant.(ParameterName) = ParameterValues(Counter);
    Parameters_DisturbedReducedClosedLoop.Plant.(ParameterName) = ParameterValues(Counter);
    
    % Simulating the Full Model
    [t_Full_1, x_Full_1] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_ClosedLoop, IC_Full, t_Disturbance, N_t, Solver);
    [t_Full_2, x_Full_2] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_DisturbedClosedLoop, x_Full_1(:,end), tf - t_Disturbance, N_t, Solver);
    t_Full = [t_Full_1, t_Disturbance + t_Full_2(2:end)];
    x_Full = [x_Full_1, x_Full_2(:,2:end)];
    x_Full_Transformed(:,:,Counter) = TransformationMatrix*x_Full;

    % Simulating the Reduced Model
    [t_Reduced_1, x_Reduced_1] = DSA(S_ReducedClosedLoop, Prop_ReducedClosedLoop, Parameters_ReducedClosedLoop, IC_Reduced, t_Disturbance, N_t, Solver);
    [t_Reduced_2, x_Reduced_2] = DSA(S_ReducedClosedLoop, Prop_ReducedClosedLoop, Parameters_DisturbedReducedClosedLoop, x_Reduced_1(:,end), tf - t_Disturbance, N_t, Solver);
    t_Reduced = [t_Reduced_1, t_Disturbance + t_Reduced_2(2:end)];
    x_Reduced(:,:,Counter) = [x_Reduced_1, x_Reduced_2(:,2:end)];   
end

%% Figure Settings
Scale = 1.5;
Figure_Width = 6.5 * Scale;
Figure_Height = 5 * Scale;
FontSize = 5 * Scale;
FontSize_Small = 4 * Scale;
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
Species = cell(Parameters_Plant.L + 3,1);
for i = 1 : Parameters_Plant.L
    if i == Parameters_Plant.L
        Species{i} = ['X_', num2str(i), '~(Regulated~Output)'];
    else
    Species{i} = ['X_', num2str(i)];
    end
end
Species{Parameters_Plant.L + 1} = ['Z^+ =~', latex(Parameters_Controller.q_tot(:,1)' * (sym('Z_', [1, Parameters_Controller.M], 'positive'))')];
Species{Parameters_Plant.L + 2} = ['Z^- =~', latex(Parameters_Controller.q_tot(:,2)' * (sym('Z_', [1, Parameters_Controller.M], 'positive'))')];
Species{Parameters_Plant.L + 3} = ['Z^0 =~', latex(Parameters_Controller.q_tot(:,3)' * (sym('Z_', [1, Parameters_Controller.M], 'positive'))')];
Handle_Tile = tiledlayout(Handle_Figure1, 2, 2);
    Handle_Tile.TileSpacing = 'compact';
    Handle_Tile.Padding = 'compact';
    xlabel(Handle_Tile, 'Time', 'FontSize', FontSize, 'FontWeight', 'bold');
    ylabel(Handle_Tile, 'Concentration', 'FontSize', FontSize, 'FontWeight', 'bold');
for i = Parameters_Plant.L : Parameters_Plant.L + 3
    Handle_Axis = nexttile;
    Handle_Axis.Box = 'on';
    Handle_Axis.BoxStyle = 'full';
    Handle_Axis.Title.Interpreter = 'latex';
    title(['$\bf{', Species{i}, '}$']);
    Handle_Axis.XLim = [0, tf];
    grid(Handle_Axis, 'on');
    hold(Handle_Axis, 'on');
    Handle_Axis.FontSize = FontSize;
    Handle_Plot = gobjects(length(ParameterValues),1);
    for Counter = 1 : length(ParameterValues)
        plot(t_Full, x_Full_Transformed(i, :, Counter), 'Color', [Colors(1,:), nthroot(Counter/length(ParameterValues), 1.5)], 'LineWidth', LineWidth_Thick);
        plot(t_Full, x_Reduced(i, :, Counter), 'Color', [Colors(2,:), nthroot(Counter/length(ParameterValues), 1.5)], 'LineWidth', LineWidth_Thick, 'LineStyle', '--');
        Handle_Plot(Counter) = plot(Handle_Axis, -t_Full, x_Full_Transformed(i, :, Counter), 'Color', [0.5, 0.5, 0.5, nthroot(Counter/length(ParameterValues), 1.5)], 'LineWidth', LineWidth);
    end
    if i == Parameters_Plant.L + 2
        Handle_Legend1 = legend(Handle_Axis, Handle_Plot, [repmat(['$\', ParameterName, ' = '], length(ParameterValues), 1), num2str(ParameterValues'), repmat('$', length(ParameterValues), 1)]);
        Handle_Legend1.Interpreter = 'latex';
        Handle_Legend1.FontSize = FontSize_Small;
        Handle_Legend1.Location = 'best';
        Handle_Legend1.AutoUpdate = 'off';
    end
    if i == Parameters_Plant.L + 3
        Handle_Legend2 = legend(Handle_Axis, 'Full Model', 'Reduced Model', 'FontSize', FontSize_Small, 'Location', 'best');
    end
end

%% Save Figure
if Save_Flag == 1
    exportgraphics(Handle_Tile, Figure1_Name, 'ContentType', 'vector', 'BackgroundColor', 'none');
end