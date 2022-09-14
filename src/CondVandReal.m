function CondVandReal()

close all

% --- OPDRACHT 5 --- %

N = 5:50;

%%% Equidistante punten in [0,1] %%%

theoCondEqPos = zeros(1, length(N));   % theoretische conditie voor de positieve equidistante punten, dus in [0,1]
condEqPos = zeros(1, length(N));       % experimentele conditie voor de positieve equidistante punten, dus in [0,1]
for k = N
    eq = linspace(0, 1, k);
    a = ones(1, k);
    for j = 1:k
        for i = 1:k
            if i ~= j
                a(j) = a(j)*(1+abs(eq(i)))/(eq(i)-eq(j));
            end
        end
    end
    theoCondEqPos(k-4) = k*max(a);
    V = Vandermonde(eq);
    condEqPos(k-4) = cond(V, 1);
end

% plotten
figure(Name="Equidistant tussen 0 en 1")
semilogy(N, condEqPos)
hold all
semilogy(N, theoCondEqPos)
% titel, labels en legende
title('Equidistant tussen 0 en 1')
xlabel('N') 
ylabel('cond(V)') 
legend('experimenteel', 'theoretisch', 'Location', 'northwest', 'FontSize', 12)
hold off

%%% Equidistante punten in [-1,1] %%%

condEq = zeros(1, length(N));          % experimentele conditie voor equidistante punten in [0,1]
for k = N
    eq = linspace(-1, 1, k);
    V = Vandermonde(eq);
    condEq(k-4) = cond(V, 1);
end

% plotten
figure(Name="Equidistant tussen -1 en 1")
semilogy(N, condEq)
% titel, labels en legende
title('Equidistant tussen -1 en 1')
xlabel('N') 
ylabel('cond(V)') 
legend('experimenteel', 'Location', 'northwest', 'FontSize', 12)

%%% Chebyshev punten %%%

condCheb = zeros(1, length(N));           % experimentele conditie voor Chebyshev punten
for k = N
    cheb = zeros(1,k);
    for v = 1:k
        cheb(v) = cos((2*v-1)*pi/(2*k));
    end
    V = Vandermonde(cheb);
    condCheb(k-4) = cond(V, 1);
end

% plotten
figure(Name="Chebyshev")
semilogy(N, condCheb)
% titel, labels en legende
title('Chebyshev')
xlabel('N') 
ylabel('cond(V)') 
legend('experimenteel', 'Location', 'northwest', 'FontSize', 12)

%%% Harmonische punten %%%

theoCondHarm = zeros(1, length(N));       % theoretische conditie voor de harmonische punten
condHarm = zeros(1, length(N));           % experimentele conditie voor de harmonische punten
for k = N
    harm = zeros(1,k);
    for v = 1:k
        harm(v) = 1/v;
    end
    a = ones(1,k);
    for j = 1:k
        for i = 1:k
            if i ~= j
                a(j) = a(j)*(1+abs(harm(i)))/(harm(i)-harm(j));
            end
        end
    end
    theoCondHarm(k-4) = k*max(a);
    V = Vandermonde(harm);
    condHarm(k-4) = cond(V, 1);
end

% plotten
figure(Name="Harmonisch")
semilogy(N,condHarm)
hold all
semilogy(N,theoCondHarm)
% titel, labels en legende
title('Harmonisch')
xlabel('N') 
ylabel('cond(V)') 
legend('experimenteel', 'theoretisch', 'Location', 'northwest', 'FontSize', 12)
hold off

%%% Vergelijking van de vier families %%%

% plotten
figure(Name="Experimentele vergelijking van de 4 families")
semilogy(N,condEqPos)
hold all
semilogy(N,condEq)
semilogy(N,condCheb)
semilogy(N,condHarm)
% titel, labels en legende
title('Vergelijking conditie van de 4 families')
xlabel('N') 
ylabel('cond(V)') 
legend('Equidistant [0,1]', 'Equidistant [-1,1]', 'Chebyshev', 'Harmonisch', 'Location', 'northwest', 'FontSize', 12)
hold off

%%% Superexponentieel gedrag bij harmonische punten %%%

% we nemen 2 punten aan het begin van de grafiek voor de conditie bij harmonische punten
x1 = N(1);
y1 = condHarm(1);
x2 = N(2);
y2 = condHarm(2);

% we stellen een stelsel op dat parameter a en b in vergelijking be^(ax) = y oplost.
% Er geldt: ln(b) + ax = ln(y) voor beide waarden van x en y
X = [x1, 1; x2, 1];
Y = [log(y1); log(y2)];
% we lossen het stelsel X*par = Y op om par = [a; ln(b)] te bekomen
par = X \ Y;
a1 = par(1);
b1 = exp(par(2));
% en maken een exponentieel verloop met de gevonden resultaten voor a en b
expHarm = zeros(1,length(N));
for k=N
    expHarm(k-N(1)+1) = b1*exp(a1*k);
end

% plotten
figure(Name="Superexponentieel gedrag bij harmonische punten")
semilogy(N,condHarm)
hold all
semilogy(N,expHarm)
% titel, labels en legende
title('Superexponentieel gedrag bij harmonische punten')
xlabel('N') 
ylabel('cond(V)') 
legend('Conditie bij harmonische punten', 'Exponentieel verloop', 'Location', 'northwest', 'FontSize', 12)
hold off

%%% (asymptotisch) exponentieel gedrag bij de overige families %%%

% [0,1]
% we nemen 2 punten aan het begin van de grafiek voor de conditie bij harmonische punten.
% Omdat we weten dat de conditie (asymptotisch) exponentieel verloopt kunnen we punten
% nemen die verder uit elkaar gelegen zijn om een accurater resultaat te bekomen
x1 = N(1);
x2 = N(10);
y1 = condEqPos(1);
y2 = condEqPos(10);

% we stellen een stelsel op dat parameter a en b in vergelijking be^(ax) = y oplost.
% Er geldt: ln(b) + ax = ln(y) voor beide waarden van x en y
X = [x1, 1; x2, 1];
Y = [log(y1); log(y2)];
% we lossen het stelsel X*par = Y op om par = [a; ln(b)] te bekomen
par = X \ Y;
a2 = par(1);
b2 = exp(par(2));
% en maken een exponentieel verloop met de gevonden resultaten voor a en b
expEqPos = zeros(1,length(N));
for k=N
    expEqPos(k-N(1)+1) = b2*exp(a2*k);
end

% plotten
figure(Name="Exponentieel gedrag bij Equidistante punten tussen 0 en 1")
semilogy(N,condEqPos)
hold all
semilogy(N,expEqPos)
% titel, labels en legende
title('Exponentieel gedrag bij Equidistante punten tussen 0 en 1')
xlabel('N') 
ylabel('cond(V)') 
legend('Conditie bij equidistante punten tussen 0 en 1', 'Exponentieel verloop', 'Location', 'northwest', 'FontSize', 12)
hold off

% [-1, 1]
% we nemen 2 punten aan het begin van de grafiek voor de conditie bij harmonische punten.
% Omdat we weten dat de conditie (asymptotisch) exponentieel verloopt kunnen we punten
% nemen die verder uit elkaar gelegen zijn om een accurater resultaat te bekomen
% de x-waarden kunnen we hergebruiken van de vorige familie punten
y1 = condEq(1);
y2 = condEq(10);

% we stellen een stelsel op dat parameter a en b in vergelijking be^(ax) = y oplost.
% Er geldt: ln(b) + ax = ln(y) voor beide waarden van x en y
X = [x1, 1; x2, 1];
Y = [log(y1); log(y2)];
% we lossen het stelsel X*par = Y op om par = [a; ln(b)] te bekomen
par = X \ Y;
a3 = par(1);
b3 = exp(par(2));
% en maken een exponentieel verloop met de gevonden resultaten voor a en b
expEq = zeros(1,length(N));
for k=N
    expEq(k-N(1)+1) = b3*exp(a3*k);
end

% plotten
figure(Name="Exponentieel gedrag bij Equidistante punten tussen -1 en 1")
semilogy(N,condEq)
hold all
semilogy(N,expEq)
% titel, labels en legende
title('Exponentieel gedrag bij Equidistante punten tussen -1 en 1')
xlabel('N') 
ylabel('cond(V)') 
legend('Conditie bij equidistante punten tussen -1 en 1', 'Exponentieel verloop', 'Location', 'northwest', 'FontSize', 12)
hold off

% Chebyshev punten
% we nemen 2 punten aan het begin van de grafiek voor de conditie bij harmonische punten.
% Omdat we weten dat de conditie (asymptotisch) exponentieel verloopt kunnen we punten
% nemen die verder uit elkaar gelegen zijn om een accurater resultaat te bekomen
% de x-waarden kunnen we hergebruiken van de vorige familie punten
y1 = condCheb(1);
y2 = condCheb(10);

% we stellen een stelsel op dat parameter a en b in vergelijking be^(ax) = y oplost.
% Er geldt: ln(b) + ax = ln(y) voor beide waarden van x en y
X = [x1, 1; x2, 1];
Y = [log(y1); log(y2)];
% we lossen het stelsel X*par = Y op om par = [a; ln(b)] te bekomen
par = X \ Y;
a4 = par(1);
b4 = exp(par(2));
% en maken een exponentieel verloop met de gevonden resultaten voor a en b
expCheb = zeros(1,length(N));
for k=N
    expCheb(k-N(1)+1) = b4*exp(a4*k);
end

% plotten
figure(Name="Exponentieel gedrag bij Chebyshev punten")
semilogy(N,condCheb)
hold all
semilogy(N,expCheb)
% titel, labels en legende
title('Exponentieel gedrag bij Chebyshev punten')
xlabel('N') 
ylabel('cond(V)') 
legend('Conditie bij Chebyshev punten', 'Exponentieel verloop', 'Location', 'northwest', 'FontSize', 12)
hold off

% zet de resultaten van de parameters in een tabel
T = table(['a'; 'b'], [a1; b1], [a2; b2], [a3; b3], [a4; b4]);
T.Properties.VariableNames = ["Parameters", "Harmonisch", "Equidistant [0, 1]", "Equidistant [-1, 1]", "Chebyshev"];
disp(T);

% --- OPDRACHT 6 --- %

%%% we kiezen de harmonische punten %%%

mean = 0;
dev = 10^(-12);
b = ones([7,1]);        % n = 7
harm = zeros(1,7);
for v = 1:7
    harm(v) = 1/v;
end
V = Vandermonde(harm);
conditie = cond(V, 1);  % conditiegetal (norm = 1)

% voorwaaartse relatieve fouten en bijhorende bovengrenzen voor tabel
relatieve_fouten = zeros(1, 5);
bovengrenzen = zeros(1, 5);

for i = 1:5
    perturbatie = normrnd(mean, dev, [1, 7]).';

    x_fault = V \ (b + perturbatie.*b);
    x_nofault = V \ (b);

    % we blijven werken met norm 1
    norm_delta = norm(perturbatie, 1);
    norm_b = norm(b, 1);

    bovengrenzen(i) = conditie*(norm_delta/norm_b);
    relatieve_fouten(i) = norm((x_fault - x_nofault),1) / norm(x_nofault, 1);
end

% tabel opstellen
T = table([1:5].', relatieve_fouten.', bovengrenzen.');
T.Properties.VariableNames = ["Perturbatie", "Relatieve fout", "Bovengrens"];
disp(T);

end