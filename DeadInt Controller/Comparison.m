%% Clear Workspace
close all
clear
clc

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 6.5 * SS;
Figure_Height = 3.3 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*0.8 * SS;
LineWidth_Thin = ScalingFactor*0.1 * SS;
MarkerSize = ScalingFactor*5 * SS;
Opacity = 0.5;
Colors = lines(10);

%% Set Figure 1
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

    
mu_1_vector = [5, 10, 15, 20];
for iter = 1 : length(mu_1_vector)
    mu_1 = mu_1_vector(iter);
%% Number of Controller Species & Charges
M = 8; 
q_tot = [   1     0     0     0     0     1     2     1; ...
            0     1     0     0     0     0     0     0; ...
            0     0     1     0     0     0     0     0; ...
            0     0     0     1     2     2     2     1; ...
        ];

%% Controller Network Stoichiometry
% Sequestration Reactions
[S_Q, Q_1, Q_2, Q_3] = SequestrationStoichiometry(1, 2, 3, M);
% Binding Reactions
[S_B, B_1, B_2, B_3] = BindingStoichiometry([4; 1; 1; 1; 8; 4], [4; 5; 6; 4; 8; 8], [5; 6; 7; 8; 7; 6], M);
% Conversion Reactions
[S_C, C_1, C_2] = ConversionStoichiometry(0, 0, M);
% Degradation Reactions
[S_D, Q_0] = DegradationStoichiometry([3; 4; 5], M);
% Actuation Reactions
h_plus = @(x_L, z, ActuationParameters) ((ActuationParameters.k*z(7) + ActuationParameters.k_prime*z(6)) / (1 + z(5)/ActuationParameters.kappa_u));
h_minus = @(x_L, z, ActuationParameters) (0);

%% Controller Network Parameters
% Sequestration Reactions
eta = 1;
% Binding Reactions
a = [10; 10; 10; 10; 10; 10] * 10;
d = [10; 10; 10; 10; 10; 10] * 10;
% Conversion Reactions
c_1 = [];
c_2 = [];
% Degradation Reactions
delta_0 = 1;
gamma = delta_0 * [1; 1; 1];
% Production Reactions;
mu = [mu_1; 0 ; 0; 5; 0; 0; 0; 0];
theta = [0; 2; 0; 0; 0; 0; 0; 0];
% Actuation Reactions
ActuationParameters.k = 1; 
ActuationParameters.k_prime = 0.5; 
ActuationParameters.kappa_u = 1;
% Dilution Reactions
delta = 0.01;

%% Select Plant
StoichiometryMatrix_Plant = StoichiometryMatrix_GeneExp();
PropensityFunction_Plant = @PropensityFunction_GeneExp;
Parameters_Plant.k_1 = 1;
Parameters_Plant.gamma_p_1 = 1;
Parameters_Plant.gamma_p_2 = 1;
L = 2;

%% Simulation Settings
tf = 100;
N_t = 1000;
Solver = 'ODE23s';
IC = zeros(L+M,1);

%% Disturbance 
DisturbanceFactor = 2;
t_Disturbance = tf/2;
DisturbedParameter = 'gamma_p_2';

%% Construct Controller
S_Controller = [zeros(M,2), eye(M), S_Q, S_C, -S_C, S_B, -S_B, S_D, -eye(M)];
PropensityFunction_Controller = @(x_L, z, h_plus, h_minus, ActuationParameters, mu, theta, eta, c_1, c_2, a, d, gamma, delta, Q_1, Q_2, C_1, C_2, B_1, B_2, B_3, Q_0) ...
                                 ([ h_plus(x_L, z, ActuationParameters); ...        % Positive Actuation
                                    h_minus(x_L, z, ActuationParameters); ...       % Negative Actuation
                                    mu + theta*x_L; ...                             % Setpoint/Sensing
                                    eta .* (Q_1*z) .* (Q_2*z); ...                  % Sequestration
                                    c_1 .* (C_1*z); ...                             % Forward Conversion         
                                    c_2 .* (C_2*z); ...                             % Backward Conversion           
                                    a .* (B_1*z) .* (B_2*z); ...                    % Forward Binding
                                    d .* (B_3*z); ...                               % Backward Binding              
                                    gamma .* (Q_0*z); ...                           % Degradation
                                    delta*z; ...                                    % Dilution
                                 ]);

%% Closed-Loop Network
S_Mutual = zeros(size(StoichiometryMatrix_Plant,1), size(S_Controller,2));
S_Mutual(1,1:2) = [1, -1];
S_ClosedLoop = [    StoichiometryMatrix_Plant,                          S_Mutual; ...
                    zeros(M, size(StoichiometryMatrix_Plant,2)),        S_Controller; ...
                ];
Parameters_ClosedLoop.Parameters_Plant = Parameters_Plant;
Parameters_ClosedLoop.L = L;
Parameters_ClosedLoop.h_plus = h_plus;
Parameters_ClosedLoop.h_minus = h_minus;
Parameters_ClosedLoop.ActuationParameters = ActuationParameters;
Parameters_ClosedLoop.mu = mu;
Parameters_ClosedLoop.theta = theta;
Parameters_ClosedLoop.eta = eta;
Parameters_ClosedLoop.c_1 = c_1;
Parameters_ClosedLoop.c_2 = c_2;
Parameters_ClosedLoop.a = a;
Parameters_ClosedLoop.d = d;
Parameters_ClosedLoop.gamma = gamma;
Parameters_ClosedLoop.delta = delta;
Parameters_ClosedLoop.Q_1 = Q_1;
Parameters_ClosedLoop.Q_2 = Q_2;
Parameters_ClosedLoop.C_1 = C_1;
Parameters_ClosedLoop.C_2 = C_2;
Parameters_ClosedLoop.B_1 = B_1;
Parameters_ClosedLoop.B_2 = B_2;
Parameters_ClosedLoop.B_3 = B_3;
Parameters_ClosedLoop.Q_0 = Q_0;
Prop_ClosedLoop = @(X, Parameters_ClosedLoop) ...
                    ([ PropensityFunction_Plant(X(1:L), Parameters_ClosedLoop.Parameters_Plant); ...
                       PropensityFunction_Controller(X(L), X(L+1:L+M), Parameters_ClosedLoop.h_plus, Parameters_ClosedLoop.h_minus, Parameters_ClosedLoop.ActuationParameters, ...
                       Parameters_ClosedLoop.mu, Parameters_ClosedLoop.theta, Parameters_ClosedLoop.eta, Parameters_ClosedLoop.c_1, Parameters_ClosedLoop.c_2, ...
                       Parameters_ClosedLoop.a, Parameters_ClosedLoop.d, Parameters_ClosedLoop.gamma, Parameters_ClosedLoop.delta, ...
                       Parameters_ClosedLoop.Q_1, Parameters_ClosedLoop.Q_2, Parameters_ClosedLoop.C_1, Parameters_ClosedLoop.C_2, Parameters_ClosedLoop.B_1, ...
                       Parameters_ClosedLoop.B_2, Parameters_ClosedLoop.B_3, Parameters_ClosedLoop.Q_0); ...
                    ]);
                
%% Disturbance
Parameters_Plant_Disturbed = Parameters_Plant;
Parameters_Plant_Disturbed.(DisturbedParameter) = Parameters_Plant.(DisturbedParameter) * DisturbanceFactor;
Parameters_ClosedLoop_Disturbed = Parameters_ClosedLoop;
Parameters_ClosedLoop_Disturbed.Parameters_Plant = Parameters_Plant_Disturbed;
                
%% Simulating the Full Model
[t_1, X_1] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_ClosedLoop, IC, t_Disturbance, N_t, Solver);
[t_2, X_2] = DSA(S_ClosedLoop, Prop_ClosedLoop, Parameters_ClosedLoop_Disturbed, X_1(:,end), tf - t_Disturbance, N_t, Solver);
t = [t_1, t_Disturbance + t_2(2:end)];
X = [X_1, X_2(:,2:end)];
z_tot = q_tot * X(L+1:L+M,:);

%% Reduced Model Construction
load ReducedModel
a_1 = a(1); a_2 = a(2); a_3 = a(3); a_4 = a(4); a_5 = a(5); a_6 = a(6);
d_1 = d(1); d_2 = d(2); d_3 = d(3); d_4 = d(4); d_5 = d(5); d_6 = d(6);
k_1 = Parameters_Plant.k_1; 
gamma_p_1 = Parameters_Plant.gamma_p_1; 
gamma_p_2 = Parameters_Plant.gamma_p_2; 
k = ActuationParameters.k;
k_prime = ActuationParameters.k_prime;
kappa_u = ActuationParameters.kappa_u;
mu_1 = mu(1); mu_4 = mu(4); theta_2 = theta(2);
F = @(t, Y, YP) f(t, Y, YP, a_1, a_2, a_3, a_4, a_5, a_6, d_1, d_2, d_3, d_4, d_5, d_6, delta, delta_0, eta, gamma_p_1, gamma_p_2, k, k_1, k_prime, kappa_u, mu_1, mu_4, theta_2);
y0est = zeros(L+M,1);
yp0est = zeros(L+M,1);
opt = odeset('RelTol', 10.0^(-7),'AbsTol',10.0^(-7));
[y0, yp0] = decic(F, 0, y0est, [], yp0est, [], opt);
[tSol1, ySol1] = ode15i(F,[0 t_Disturbance], y0, yp0, opt);
gamma_p_2 = gamma_p_2 * DisturbanceFactor;
F = @(t, Y, YP) f(t, Y, YP, a_1, a_2, a_3, a_4, a_5, a_6, d_1, d_2, d_3, d_4, d_5, d_6, delta, delta_0, eta, gamma_p_1, gamma_p_2, k, k_1, k_prime, kappa_u, mu_1, mu_4, theta_2);
y0est = ySol1(end,:)';
[y0, yp0] = decic(F, 0, y0est, [], yp0est, [], opt);
[tSol2, ySol2] = ode15i(F,[t_Disturbance, tf], y0, yp0, opt);
tSol = [tSol1; tSol2(2:end)];
ySol = [ySol1; ySol2(2:end,:)];
    
%% Axis for Simulations
% set (0, 'CurrentFigure' , Handle_Figure1)
% Handle_Axis1 = gobjects(6,1);
% for i = 1 : 6
%     Handle_Axis1(i) = subplot(3,2,i);
%         Handle_Axis1(i).Box = 'on';
%         Handle_Axis1(i).BoxStyle = 'full';
%         Handle_Axis1(i).LineWidth = LineWidth_Thin;
%         Handle_Axis1(i).FontSize = FontSize;
%         hold(Handle_Axis1(i), 'on');
%         grid(Handle_Axis1(i), 'on');
%         Handle_Axis1(i).XMinorGrid = 'on';
%         Handle_Axis1(i).YMinorGrid = 'on';
%         switch i
%             case 1
%                 Handle_Axis1(i).Title.String = ['$x_1$', '(Input)'];
%                 Handle_Axis1(i).XTickLabel = [];
% %                 Handle_Axis1(i).Position = [0.05, 0.56, 0.86, 0.38];
%             case 2
%                 Handle_Axis1(i).Title.String = ['$x_2$', '(Output)'];
%                 Handle_Axis1(i).XTickLabel = [];
% %                 Handle_Axis1(i).Position = [0.56, 0.56, 0.43, 0.38];
%             case 3
%                 Handle_Axis1(i).Title.String = '$z^+$';
%                 Handle_Axis1(i).XTickLabel = [];
% %                 Handle_Axis1(i).Position = [0.56, 0.56, 0.43, 0.38];
%             case 4
%                 Handle_Axis1(i).Title.String = '$z^-$';
%                 Handle_Axis1(i).XLabel.String = [];
% %                 Handle_Axis1(i).Position = [0.05, 0.1, 0.43, 0.38];
%             case 5
%                 Handle_Axis1(i).Title.String = '$z^0$';
%                 Handle_Axis1(i).XLabel.String = 'Time';
% %                 Handle_Axis1(i).Position = [0.05, 0.1, 0.43, 0.38];
%             case 6
%                 Handle_Axis1(i).Title.String = '$z^\star$';
%                 Handle_Axis1(i).XLabel.String = 'Time';
% %                 Handle_Axis1(i).Position = [0.56, 0.1, 0.43, 0.38];
%         end
%         Handle_Axis1(i).Title.Interpreter = 'latex';
%         Handle_Axis1(i).XLim = [0, tf];
%         if i <= 2
%             plot(Handle_Axis1(i), t, X(i, :), 'Color', Colors(1,:), 'LineWidth', LineWidth_Thick);
%             plot(Handle_Axis1(i), tSol, ySol(:,i), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick, 'LineStyle', '--');
%         else
%             plot(Handle_Axis1(i), t, z_tot(i-2,:), 'Color', Colors(1,:), 'LineWidth', LineWidth_Thick);
%             plot(Handle_Axis1(i), tSol, ySol(:,i), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick, 'LineStyle', '--');
%         end
%             
% end

i = 1;
    if iter == 1
        Handle_Axis1(i) = subplot(1,1,i);
        box(Handle_Axis1(i), 'on');
    end
        hold(Handle_Axis1(i), 'on');
        Handle_Axis1(i).Position = [0.1, 0.1, 0.87, 0.87];
        Handle_Axis1(i).Box = 'on';
        Handle_Axis1(i).BoxStyle = 'full';
        Handle_Axis1(i).LineWidth = LineWidth_Thin;
        Handle_Axis1(i).FontSize = FontSize;
        hold(Handle_Axis1(i), 'on');
        grid(Handle_Axis1(i), 'on');
        Handle_Axis1(i).XMinorGrid = 'on';
        Handle_Axis1(i).YMinorGrid = 'on';
        Handle_Axis1(i).YLim = [0, 18];
        switch i
            case 1
%                 Handle_Axis1(i).XTickLabel = [];
                Handle_Axis1(i).XLabel.String = 'Time';
                Handle_Axis1(i).YLabel.String = 'Regulated Output';
%                 Handle_Axis1(i).Position = [0.56, 0.56, 0.43, 0.38];
        end
        Handle_Axis1(i).XLim = [0, tf];
        Handle_Axis1(i).XLabel.Position(2) = 1.6;
        h_1 = plot(Handle_Axis1(i), t, X(2, :), 'Color', [Colors(1,:), mu(1)/20], 'LineWidth', LineWidth_Thick);
     	h_2 = plot(Handle_Axis1(i), tSol, ySol(:,2), 'Color', [Colors(2,:), mu(1)/20], 'LineWidth', LineWidth_Thick, 'LineStyle', '--');
        Handle_Axis1(i).XRuler.Axle.LineWidth = 3;
        Handle_Axis1(i).YRuler.Axle.LineWidth = 3;
end
legend(Handle_Axis1(1), [h_1, h_2], {'Full Model', 'Reduced Model'});
annotation(Handle_Figure1,'textbox',...
    [0.799041666666667 0.528754629068906 0.194135423554315 0.11406481537554],...
    'String',{'$\mu_1 = 20$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Helvetica Neue');

annotation(Handle_Figure1,'textbox',...
    [0.799041666666667 0.417041666105942 0.194135423554314 0.11406481537554],...
    'String',{'$\mu_1 = 15$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(Handle_Figure1,'textbox',...
    [0.799041666666667 0.302388888328166 0.194135423554314 0.11406481537554],...
    'String',{'$\mu_1 = 10$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(Handle_Figure1,'textbox',...
    [0.799041666666667 0.187736110550388 0.172968689600626 0.11406481537554],...
    'String',{'$\mu_1 = 5$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');



%% Save Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'DeadInts_Simulations', '-dpdf', '-painters');
end


