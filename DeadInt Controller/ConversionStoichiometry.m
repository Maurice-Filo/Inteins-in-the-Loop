function [S_C, C_1, C_2] = ConversionStoichiometry(Species_1, Species_2, M)
    if Species_1 == 0
        S_C = zeros(0,M);
        C_1 = zeros(0,M);
        C_2 = zeros(0,M);
    else
        p_0 = size(Species_1,1);
        C_1 = zeros(p_0,M);
        C_2 = zeros(p_0,M);
        keyboard
        for l = 1 : p_0
            C_1(l, Species_1(l)) = 1;
            C_2(l, Species_2(l)) = 1;
        end
        S_C = (C_2 - C_1)';
    end
end

