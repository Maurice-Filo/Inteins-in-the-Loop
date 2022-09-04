function [S_D, Q_0] = Degradation(Rxns, M)
%% Extract Species
m_0 = length(Rxns);
Species = zeros(m_0,1);
for i = 1 : length(Rxns)
    Rxn = Rxns{i};
    Rxn = Rxn(~isspace(Rxn));
    Index = strfind(Rxn, '->');
    Reactant = Rxn(1:Index-1);
 	Reactant([strfind(Reactant, 'Z'), strfind(Reactant, '_')]) = [];
	Species(i) = str2double(Reactant);
end

%% Construct Matrices
    if isempty(Species)
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

