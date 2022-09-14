function A = cubicsplinematrix(t, tbefore, tafter, periodic)
% t is een vector met knooppunten t_0 ... t_n. Deze zijn geordend en allen
% verschillend
% tbefore is een (geordende) vector met extra knooppunten die we vooraan toevoegen
% tafter is een (geordende) vector met extra knooppunten die we achteraan toevoegen
% periodic is een logische constante. Een periodic = 1 komt overeen met periodische
% voorwaarden. Een periodic = 0 komt overeen met natuurlijke voorwaarden.
% A is de matrix die de interpolatievoorwaarden beschrijft.


n = length(t) - 1;
A = zeros(n+3);
v = [tbefore, t, tafter];
orde = 4;      % want graad is 3

%%% bereken 2de tot (n+2)de rij van A %%%

for i = 2:n+2
    for j = (i-1):(i+1)
        pp = bspline(v(j:j+4));
        A(i, j) = ppval(pp, v(i+2));
    end
end

%%% randvoorwaarden (eerste en laatste rij van A) %%%

% eerte rij
for j = 1:3         
        pp = bspline(v(j:j+4));
        dp = fnder(pp, 2);
        A(1, j) = ppval(dp, v(orde));
end

if periodic == 0    % natuurlijke voorwaarden
    % laatse rij
    for j = n+1:n+3
        pp = bspline(v(j:j+4));
        dp = fnder(pp, 2);
        A(n+3, j) = ppval(dp, v(orde+n));
    end
else                % periodische voorwaarden
    % eerste rij 2de deel
    for j = n+1:n+3
        pp = bspline(v(j:j+4));
        dp = fnder(pp, 2);
        A(1, j) = -ppval(dp, v(orde+n));
    end
    % laatste rij 1ste deel
    for j = 1:3
        pp = bspline(v(j:j+4));
        dp = fnder(pp, 1);
        A(n+3, j) = ppval(dp, v(orde));
    end
    % laatste rij 2de deel
    for j = n+1:n+3
        pp = bspline(v(j:j+4));
        dp = fnder(pp, 1);
        A(n+3, j) = -ppval(dp, v(orde+n));
    end
end

end