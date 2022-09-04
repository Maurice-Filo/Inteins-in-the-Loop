function [time_vector, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver)
time_vector = linspace(0, tf, N);
switch Solver
    case 'ODE45'
        Options = odeset('NonNegative',1, 'RelTol', 1e-5);
        [~, x] = ode45(@(t, x) RHS(t, x, StoichiometryMatrix, PropensityFunction, NetworkParameters), time_vector, x0, Options);
    case 'ODE15s'
        Options = odeset('NonNegative',1, 'RelTol', 1e-5);
        [~, x] = ode15s(@(t, x) RHS(t, x, StoichiometryMatrix, PropensityFunction, NetworkParameters), time_vector, x0, Options);
    case 'ODE23s'
        Options = odeset('AbsTol', 1e-5, 'RelTol', 1e-5);
        [~, x] = ode23s(@(t, x) RHS(t, x, StoichiometryMatrix, PropensityFunction, NetworkParameters), time_vector, x0, Options);
end
x = x';
end

function dx = RHS(~, x, StoichiometryMatrix, PropensityFunction, NetworkParameters)
    dx = StoichiometryMatrix * PropensityFunction(x, NetworkParameters);
end