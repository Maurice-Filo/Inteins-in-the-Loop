function [S_Q, Q_1, Q_2, Q_3] = SequestrationStoichiometry(Species_1, Species_2, Species_3, M)
    if Species_1 == 0
        S_Q = [];
        Q_1 = [];
        Q_2 = [];
        Q_3 = [];
    else
        m = size(Species_1,1);
        Q_1 = zeros(m,M);
        Q_2 = zeros(m,M);
        Q_3 = zeros(m,M);
        for l = 1 : m
            Q_1(l, Species_1(l)) = 1;
            Q_2(l, Species_2(l)) = 1;
            for j = 1 : size(Species_3,2)
                if Species_3(l,j) ~= 0
                    Q_3(l, Species_3(l,j)) = 1;
                end
            end
        end
        S_Q = (Q_3 - Q_2 - Q_1)';
    end
end

