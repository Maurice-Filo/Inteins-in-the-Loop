function [S_C, C_1, C_2] = Conversion(Rxns, M)
%% Extract Species
p_0 = length(Rxns);
Species_1 = zeros(p_0,1);
Species_2 = zeros(p_0,1);
for i = 1 : length(Rxns)
    Rxn = Rxns{i};
    Rxn = Rxn(~isspace(Rxn));
    Index = strfind(Rxn, '<->');
    Reactant = Rxn(1:Index-1);
    Product = Rxn(Index+3:end);
 	Reactant([strfind(Reactant, 'Z'), strfind(Reactant, '_')]) = [];
	Species_1(i) = str2double(Reactant);
    Product([strfind(Product, 'Z'), strfind(Product, '_')]) = [];
    Species_2(i) = str2double(Product);
end

%% Construct Matrices
    if isempty(Species_1)
        S_C = zeros(0,M);
        C_1 = zeros(0,M);
        C_2 = zeros(0,M);
    else
        p_0 = size(Species_1,1);
        C_1 = zeros(p_0,M);
        C_2 = zeros(p_0,M);
        for l = 1 : p_0
            C_1(l, Species_1(l)) = 1;
            C_2(l, Species_2(l)) = 1;
        end
        S_C = (C_2 - C_1)';
    end
end

