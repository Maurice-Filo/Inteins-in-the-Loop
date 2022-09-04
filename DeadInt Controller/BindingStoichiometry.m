function [S_B, B_1, B_2, B_3] = BindingStoichiometry(Species_1, Species_2, Species_3, M)
    if Species_1 == 0
        S_B = [];
        B_1 = [];
        B_2 = [];
        B_3 = [];
    else
        p = size(Species_1,1);
        B_1 = zeros(p,M);
        B_2 = zeros(p,M);
        B_3 = zeros(p,M);
        for l = 1 : p
            B_1(l, Species_1(l)) = 1;
            B_2(l, Species_2(l)) = 1;
            B_3(l, Species_3(l)) = 1;
        end
        S_B = (B_3 - B_2 - B_1)';
    end
end

