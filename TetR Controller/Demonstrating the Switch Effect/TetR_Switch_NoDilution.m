%% Clear Workspace
close all
clear
clc

k_prime_vector = [0, 0.5, 1, 2, 3];
for iter = 1 : length(k_prime_vector)
%% Select Plant: {'Birth-Death', 'Gene Expression', 'Six Node'}
PlantChoice = 'Gene Expression';
switch PlantChoice
    % Birth-Death
    case 'Birth-Death'
        S_Plant = StoichiometryMatrix_BD();
        Prop_Plant = @PropensityFunction_BD;
        Parameters_Plant.gamma_1 = 0.1;
        Parameters_Plant.L = 1;
        Figure1_Name = 'Switch_TetR_BD.pdf';
        % Swept Parameter
        ParameterName = 'gamma_1';
        ParameterValues = [0.1, 0.5, 1, 1.5];
        
   	% Gene Expression
    case 'Gene Expression'
        S_Plant = StoichiometryMatrix_GeneExp();
        Prop_Plant = @PropensityFunction_GeneExp;
        Parameters_Plant.gamma_1 = 1;
        Parameters_Plant.gamma_2 = 1;
        Parameters_Plant.k_1 = 1;
        Parameters_Plant.L = 2;
        Figure1_Name = 'Switch_TetR_GeneExp.pdf';
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
        Figure1_Name = 'Switch_TetR_SixNode.pdf';
        % Swept Parameter
        ParameterName = 'gamma_6';
        ParameterValues = [0.1, 0.5, 1, 1.5];
end

%% Select Controller
[Parameters_Controller, S_Controller, Prop_Controller] = Design_TetRController();
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
    Parameters_Controller.Parameters_Actuation.k = 0.1;
    Parameters_Controller.Parameters_Actuation.k_prime = k_prime_vector(iter);
    Parameters_Controller.Parameters_Actuation.kappa_u = 1;
S_ReducedController = StoichiometryMatrix_ReducedController();
Prop_ReducedController = @PropensityFunction_ReducedController;

%% Select Reduced Control Action
Parameters_ReducedController.ControlAction = @ReducedControlAction_TetRController;

%% Simulation Settings
tf = 500;
N_t = 1000;
Solver = 'ODE23s';

%% Initial Conditions
IC_Full = zeros(Parameters_Controller.M + Parameters_Plant.L, 1);
TransformationMatrix = [eye(Parameters_Plant.L), zeros(Parameters_Plant.L, Parameters_Controller.M); ...
                        zeros(3, Parameters_Plant.L), Parameters_Controller.q_tot'];
IC_Reduced = TransformationMatrix*IC_Full;

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

%% Simulating the Full Model
[t_Full, x_Full] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_ClosedLoop, IC_Full, tf, N_t, Solver);
x_Full_Transformed = TransformationMatrix*x_Full;

%% Simulating the Reduced Model
[t_Reduced, x_Reduced] = DSA(S_ReducedClosedLoop, Prop_ReducedClosedLoop, Parameters_ReducedClosedLoop, IC_Reduced, tf, N_t, Solver);

%% Extract Control Action
u = zeros(1, length(t_Reduced));
for j = 1 : length(t_Reduced)
    [u(j), ~] = Parameters_ReducedController.ControlAction(0, x_Reduced(Parameters_Plant.L+1, j), 0, x_Reduced(Parameters_Plant.L+3, j), Parameters_ReducedController); 
end

%% Compute Control Map
z_P_vector = linspace(0, max(x_Reduced(Parameters_Plant.L+1,:)), 50);
z_0_vector = linspace(0, max(x_Reduced(Parameters_Plant.L+3,:)), 50);
ControlMap = zeros(length(z_P_vector), length(z_0_vector));
for j = 1 : length(z_P_vector)
    for k = 1 : length(z_0_vector)
        [ControlMap(j,k), ~] = Parameters_ReducedController.ControlAction(0, z_P_vector(j), 0, z_0_vector(k), Parameters_ReducedController); 
    end
end

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure1_Width = 2.65 * SS;
Figure1_Height = 1.25 * SS;
Figure2_Width = 1.5 * SS;
Figure2_Height = 1.1 * SS;
FontSize = ScalingFactor*3 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.45 * SS;
LineWidth_Thick = ScalingFactor*1 * SS;
LineWidth_Thin = ScalingFactor*0.1 * SS;
MarkerSize = ScalingFactor*5 * SS;
Opacity = 0.5;
Colors = lines(10);
    
%% Set Figure 1
if iter == 1
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure1_Width, Figure1_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Axis for Simulations
set (0, 'CurrentFigure' , Handle_Figure1)
Handle_Axis1 = gobjects(4,1);
i = 1;
    Handle_Axis1(i) = subplot(2,2,i);
        Handle_Axis1(i).Box = 'on';
        Handle_Axis1(i).BoxStyle = 'full';
        Handle_Axis1(i).LineWidth = LineWidth_Thin;
        Handle_Axis1(i).FontSize = FontSize;
        hold(Handle_Axis1(i), 'on');
        grid(Handle_Axis1(i), 'on');
        Handle_Axis1(i).XMinorGrid = 'off';
        Handle_Axis1(i).YMinorGrid = 'off';
        switch i
            case 1
                Handle_Axis1(i).Title.String = 'Regulated Output';
                Handle_Axis1(i).YLabel.String = '$x_2$';
                Handle_Axis1(i).YLabel.Interpreter = 'latex';
                Handle_Axis1(i).XLabel.String = 'Time';
                Handle_Axis1(i).YLim = [0, 11];
                Handle_Axis1(i).Position = [0.12, 0.24, 0.84, 0.64];
            case 2
                Handle_Axis1(i).Title.String = '$z^+ = z_1 + 2z_5 + z_6$';
                Handle_Axis1(i).XTickLabel = [];
                Handle_Axis1(i).Position = [0.56, 0.56, 0.43, 0.38];
            case 3
                Handle_Axis1(i).Title.String = '$z^- = z_2$';
                Handle_Axis1(i).XLabel.String = 'Time';
                Handle_Axis1(i).Position = [0.05, 0.1, 0.43, 0.38];
            case 4
                Handle_Axis1(i).Title.String = '$z^0 = z_3 + 2z_4 + z_6$';
                Handle_Axis1(i).XLabel.String = 'Time';
                Handle_Axis1(i).Position = [0.56, 0.1, 0.43, 0.38];
        end
        Handle_Axis1(i).Title.Interpreter = 'latex';
        Handle_Axis1(i).XLim = [0, tf];
annotation(Handle_Figure1,'textarrow',[0.328621908127208 0.162544169611307],...
    [0.612676056338028 0.838028169014085],...
    'Color',[0.301960784313725 0.745098039215686 0.933333333333333],...
    'LineWidth',8,...
    'Interpreter','latex',...
    'HeadWidth',20,...
    'HeadLength',20,...
    'FontSize',15);
annotation(Handle_Figure1,'textbox',...
    [0.330128205128205 0.445478873911038 0.170841174247937 0.222535210596004],...
    'String',{'$k'' \uparrow$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue');

end
plot(Handle_Axis1(i), t_Full, x_Full_Transformed(Parameters_Plant.L + i - 1, :), 'Color', [0, 0, 0, iter/length(k_prime_vector)], 'LineWidth', LineWidth);

%% Figures for Control Maps
if iter == 1
    Handle_Figure2 = gobjects(length(k_prime_vector),1);
end
Handle_Figure2(iter) = figure();
    Handle_Figure2(iter).Color = [1 1 1];
    Handle_Figure2(iter).PaperUnits = 'centimeters';
    Handle_Figure2(iter).Units = 'centimeters';
    Handle_Figure2(iter).Position = [0, 0, Figure2_Width, Figure2_Height];
    Handle_Figure2(iter).PaperPositionMode = 'auto';
    Handle_Figure2(iter).PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];  
    
%% Axes for Control Maps
if iter == 1
    Handle_Axis2 = gobjects(length(k_prime_vector),1);
end
Handle_Axis2(iter) = axes(Handle_Figure2(iter)); 
    Handle_Axis2(iter).Position = [0.23, 0.21, 0.56, 0.78];
    Handle_Axis2(iter).Box = 'on';
    Handle_Axis2(iter).BoxStyle = 'full';
    Handle_Axis2(iter).LineWidth = LineWidth_Thin;
    Handle_Axis2(iter).FontSize = FontSize;
    Handle_Axis2(iter).XLim = [z_P_vector(1), z_P_vector(end)];
    Handle_Axis2(iter).YLim = [z_0_vector(1), z_0_vector(end)];
%     Handle_Axis2(iter).XTick = [z_plus_vector(1), round(z_plus_vector(end), 1, 'Significant')];
%     Handle_Axis2(iter).YTick = [z_0_vector(1), round(z_0_vector(end), 1, 'Significant')];
    Handle_Axis2(iter).ZTick = [0, 10, 20, 30, 40, 50];
    hold(Handle_Axis2(iter), 'on');
    grid(Handle_Axis2(iter), 'on');
    Handle_Axis2(iter).XMinorGrid = 'off';
    Handle_Axis2(iter).YMinorGrid = 'off';
    Handle_Axis2(iter).XLabel.String = '$z^+$';
    Handle_Axis2(iter).YLabel.String = '$z^0$';
    Handle_Axis2(iter).ZLabel.String = '$u = \mathcal U(z^+, z^0)$';
    Handle_Axis2(iter).XLabel.Interpreter = 'latex';
    Handle_Axis2(iter).YLabel.Interpreter = 'latex';
    Handle_Axis2(iter).ZLabel.Interpreter = 'latex';
    Handle_Axis2(iter).ZScale = 'linear';
    Handle_Axis2(iter).ColorScale = 'linear';
    Handle_Axis2(iter).View = [-70, 20];
%     Handle_Axis2(iter).XLabel.Position = [41.9227  227.0416    0.0001];
%     Handle_Axis2(iter).YLabel.Position = [-6.1223  497.1394    0.0013];
surf(Handle_Axis2(iter), z_P_vector, z_0_vector, ControlMap', 'FaceAlpha', 0.8);
shading interp
colormap(Handle_Axis2(iter), turbo);
plot3(Handle_Axis2(iter), x_Reduced(Parameters_Plant.L+1,:), x_Reduced(Parameters_Plant.L+3,:), u, 'Color', [0, 0, 0], 'LineWidth', LineWidth_Thick);

[XX, YY] = meshgrid(z_P_vector, z_0_vector);
ZZ = ControlMap';
spacing = 5; 
for k = 1 : spacing : length(XX(:,1))
    plot3(XX(:,k), YY(:,k), ZZ(:,k), '-k', 'LineWidth', 0.1);
    plot3(XX(k,:), YY(k,:), ZZ(k,:), '-k', 'LineWidth', 0.1);
end
end

Handle_Legend = legend(Handle_Axis1(1), {['$k'' = ', num2str(k_prime_vector(1)), '$'], ...
                                         ['$k'' = ', num2str(k_prime_vector(2)), '$'], ...
                                         ['$k'' = ', num2str(k_prime_vector(3)), '$'], ...
                                         ['$k'' = ', num2str(k_prime_vector(4)), '$'], ...
                                         ['$k'' = ', num2str(k_prime_vector(5)), '$']});
Handle_Legend.Interpreter = 'latex';
Handle_Legend.Location = 'best';

%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'SwitchingSimulation', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
end

if Save_Flag == 1  
    for iter = 1 : length(k_prime_vector)
        print(Handle_Figure2(iter), ['ControlMap_', num2str(iter)], '-dpng', '-r300');
    end
end
