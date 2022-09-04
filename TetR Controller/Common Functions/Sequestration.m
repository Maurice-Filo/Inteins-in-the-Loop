function [S_Q, Q_1, Q_2, Q_3] = Sequestration(Rxns, M)
%% Extract Species
m = length(Rxns);
Species_12 = zeros(m,2);
Species_3 = zeros(m,10);
k_max = 1;
for i = 1 : length(Rxns)
    Rxn = Rxns{i};
    Rxn = Rxn(~isspace(Rxn));
    Index = strfind(Rxn, '->');
    Reactants = split(Rxn(1:Index-1), '+');
    Products = split(Rxn(Index+2:end), '+');
    for k = 1 : 2
        Reactant = Reactants{k};
        Reactant([strfind(Reactant, 'Z'), strfind(Reactant, '_')]) = [];
        Species_12(i,k) = str2double(Reactant);
    end
    if length(Products) > k_max
        k_max = length(Products);
    end
    for k = 1 : length(Products)
        Product = Products{k};
        Product([strfind(Product, 'Z'), strfind(Product, '_')]) = [];
        Species_3(i,k) = str2double(Product);
    end
end
Species_1 = Species_12(:,1);
Species_2 = Species_12(:,2);
Species_3(:,k_max+1:end) = [];

%% Construct Matrices
    if isempty(Species_1)
        S_Q = [];
        Q_1 = [];
        Q_2 = [];
        Q_3 = [];
    else
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

