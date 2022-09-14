function [c,y] = cubicsplinesolve(t, tbefore, tafter, periodic, f, i)
% t is een vector met knooppunten t_0 ... t_n. Deze zijn geordend en allen
% verschillend
% tbefore is een (geordende) vector met extra knooppunten die we vooraan toevoegen
% tafter is een (geordende) vector met extra knooppunten die we achteraan toevoegen
% periodic is een logische constante. Een periodic = 1 komt overeen met periodische
% voorwaarden. Een periodic = 0 komt overeen met natuurlijke voorwaarden.
% f is een vector met de interpolatiewaarden geassocieerd aan de knooppunten in t
% c is de oplossing van het stelsel Ac = f. Het geeft de coëfficiënten van de
% interpolerende spline s ten opzichte van de b-spline-basis
% y is een vector die de gevonden interpolerende spline s evalueert in k
% equidistante punten tussen t_0 en t_n.

A = cubicsplinematrix(t, tbefore, tafter, periodic);
c = A \ f;

x = linspace(t(1), t(length(t)), i);
t = [tbefore, t, tafter];
% extra waarden voor en na t om t in de randwaarden juist te kunnen evalueren
tExtra = [t(1)-3; t(1)-2; t(1)-1; t(:); t(length(t))+1; t(length(t))+2; t(length(t))+3];
% extra nullen als coëfficiënten indien punten buiten t geëvalueerd moeten worden
cExtra = [0; 0; 0; c(:); 0; 0; 0];

y = zeros(1, i);
for i = 1:length(x)
    % j is de eerste index voor t waarvoor geldt: x(i) >= t(j)
    j = find(x(i) >= t,1,'last');
    if isempty(j)
        % indien j niet gevonden, dan geldt: x(i) < t(1)
        % en stellen we j gelijk aan 0 om hierop te duiden
        j = 0;
    end
    if j == 0 || j == length(t)
        % de te evalueren waarde ligt buiten t
        % dus kunnen we y(i) gelijk stellen aan 0
        y(i) = 0;
    else
        % de te evalueren waarde ligt binnen de grenzen van t
        % de coëfficiënten die nodig zijn om x(i) te evalueren:
        cEval = cExtra(j:j+3);
        % de t-waarden die nodig zijn om x(i) te evalueren:
        tEval = tExtra(j:j+6);
        for l = 1:3
            % op basis van de oefenzitting:
            % c[l]_m voor m=j-3+i...j
            for m = 4:-1:l+1
                % recursiebetrekking (voor B-splines) in dalende volgorde,
                % zodat we enkel die c's overschrijven die in de l-de stap
                % niet meer gebruikt worden.
                a = (x(i)-tEval(m))/(tEval(m+4-l)-tEval(m));    % berekenen van de alpha waarde
                cEval(m) = a*cEval(m) + (1-a)*cEval(m-1);       % recursiebetrekking zelf
            end
        end
        % het resultaat voor y(i) staat na de recursiebetrekking in cEval(4)
        y(i) = cEval(4);
    end
end

end