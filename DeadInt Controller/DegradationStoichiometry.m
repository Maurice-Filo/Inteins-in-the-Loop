function [S_D, Q_0] = DegradationStoichiometry(Species, M)
    if Species == 0
        S_D = [];
        Q_0 = [];
    else
        m_0 = size(Species,1);
        Q_0 = zeros(m_0,M);
        for l = 1 : m_0
            Q_0(l, Species(l)) = 1;
        end
        S_D = -Q_0';
    end
end

