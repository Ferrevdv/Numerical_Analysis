function CondVandImag()

close all

% --- OPDRACHT 7 --- %

N = 1:50;

%%% Conditie voor n-de eenheidswortels (niet-inductief) %%%

con = zeros(1, length(N));           % conditie voor N-de eenheidswortels
for k = N
    points = zeros(1, k);            % k-de eenheidswortels
    for v = 1:k
        points(v) = exp(2*pi*1i*(v-1)/k);
    end
    V = Vandermonde(points);
    con(k) = cond(V,2);
end

% plotten
figure(Name="Conditie voor n-de eenheidswortels (niet-inductief)")
plot(N,con)
% titel en labels
title('Conditie voor n-de eenheidswortels (niet-inductief)')
xlabel('N')
ylabel('cond(V)')

% --- OPDRACHT 8 --- %

%%% 6de eenheidswortels %%%

k = 6;
points6 = zeros(1, k);               % 6de eenheidswortels
for v = 1:k
    points6(v) = exp(2*pi*1i*(v-1)/k);
end

%%% 7de eenheidswortels %%%

k = 7;
points7 = zeros(1, k);               % 7de eenheidswortels
for v = 1:k
    points7(v) = exp(2*pi*1i*(v-1)/k);
end

% plotten
figure(Name="Vergelijking 6de en 7de eenheidswortels")
axis on
plot(points6, 'g*', 'Color', 'b')
hold all
plot(points7, 'bo', 'Color', 'r')
xline(0)    % teken de reële as
yline(0)    % teken de imaginaire as
% titel, labels en legende
title('Vergelijking 6de en 7de eenheidswortels')
xlabel('Re')
ylabel('Im')
legend('6de eenheidswortels', '7de eenheidswortels', 'Location', 'northwest', 'FontSize', 12)
hold off

% --- OPDRACHT 9 --- %

%%% 8ste quasi-cyclische eenheidswortels (inductief) %%%

n = 8;
k = 0;
pointsCycl = ones(1, n);             % n-de quasi-cyclische eenheidswortels (met n = 8)
for v = 2:n     % 1ste punt al geïnitialiseerd op 1
    while 2^(k+1) < v
        k = k + 1;
    end
    pointsCycl(v) = exp(2*pi*1i*(2*(v-2^k)-1)/(2^(k+1)));
end

% plotten
figure(Name="8ste quasi-cyclische eenheidswortels (inductief)")
plot(pointsCycl, 'bo')
xline(0)    % teken de reële as
yline(0)    % teken de imaginaire as
% titel en labels
title('8ste quasi-cyclische eenheidswortels (inductief)')
xlabel('Re')
ylabel('Im')

% --- OPDRACHT 10 --- %

%%% Conditie voor de n-de quasi-cyclische eenheidswortels %%%

N = 1:64;

conQC = zeros(1, length(N));         % Conditie voor de n-de quasi-cyclische eenheidswortels
for n = N
    k = 0;
    pointsCycl = ones(1, n);         % n-de quasi-cyclische eenheidswortels
    for v = 2:n     % 1ste wortel al geïnitialiseerd op 1
        while 2^(k+1) < v
            k = k + 1;
        end
        pointsCycl(v) = exp(2*pi*1i*(2*(v-2^k)-1)/(2^(k+1)));
    end
    V = Vandermonde(pointsCycl);
    conQC(n) = cond(V, 2);
end

% plotten
figure(Name="Conditie n-de quasi-cyclische eenheidswortels (inductief)")
semilogy(N,conQC)
% titel en labels
title('Conditie n-de quasi-cyclische eenheidswortels')
xlabel('N')
ylabel('cond(V)')

%%% Conditie voor Van der Corput knooppunten %%%

N = 1:64;
conVdC = zeros(1,length(N));        % Conditie voor Van der Corput knooppunten
for k = N
    pointsVdC = ones(1,n);          % eerste k Van der Corput knooppunten
    for n = 1:k-1
        bin = dec2bin(n);
        ntilde = 0;
        for j = 1:length(bin)
            ntilde = ntilde + str2double(bin(length(bin)+1-j))*2^(-j);
        end
        pointsVdC(n+1) = exp(2*pi*1i*ntilde);
    end
    V = Vandermonde(pointsVdC);
    conVdC(k) = cond(V,2);
end

% plotten
figure(Name="Conditie Van der Corput knooppunten")
plot(N,conVdC)
% titel en labels
title('Conditie Van der Corput knooppunten')
xlabel('N')
ylabel('cond(V)')

end