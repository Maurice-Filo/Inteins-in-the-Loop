function [S_B, B_1, B_2, B_3] = Binding(Rxns, M)
%% Extract Species
p = length(Rxns);
Species_12 = zeros(p,2);
Species_3 = zeros(p,1);
for i = 1 : length(Rxns)
    Rxn = Rxns{i};
    Rxn = Rxn(~isspace(Rxn));
    Index = strfind(Rxn, '<->');
    Reactants = split(Rxn(1:Index-1), '+');
    Products = split(Rxn(Index+3:end), '+');
    for k = 1 : 2
        Reactant = Reactants{k};
        Reactant([strfind(Reactant, 'Z'), strfind(Reactant, '_')]) = [];
        Species_12(i,k) = str2double(Reactant);
    end
    k = 1;
    Product = Products{k};
    Product([strfind(Product, 'Z'), strfind(Product, '_')]) = [];
	Species_3(i,k) = str2double(Product);
end
Species_1 = Species_12(:,1);
Species_2 = Species_12(:,2);

%% Construct Matrices
    if isempty(Species_1)
        B_1 = zeros(0,M);
        B_2 = zeros(0,M);
        B_3 = zeros(0,M);
        S_B = (B_3 - B_2 - B_1)';
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

