close all

% --- OPDRACHT 15 --- %

% lissajous waarden voor 200 punten
[gammax_1, gammay_1] = lissajous(2, 3, pi/4, 200);
[gammax_2, gammay_2] = lissajous(3/2, 5/2, pi/2, 200);

l=20;       % l interpolatiepunten
k=200;      % spline geëvalueerd in k punten

% t, tbefore en tafter
t = linspace(0, 2*pi, l);
tbefore = [-3, -2, -1];
tafter = [7, 8, 9];

% l interpolatiepunten selecteren
[lisx_1, lisy_1] = lissajous(2, 3, pi/4, l);
[lisx_2, lisy_2] = lissajous(3/2, 5/2, pi/2, l);

for i = 0:1
    % volgend proces 2 keer herhalen, 1 keer voor elk van de
    % lissajousfiguren, met de bijhorende parameters
    if i == 0
        gammax = gammax_1;
        gammay = gammay_1;
        lisx = lisx_1;
        lisy = lisy_1;
    else
        gammax = gammax_2;
        gammay = gammay_2;
        lisx = lisx_2;
        lisy = lisy_2;
    end

    % lissajousfiguren plotten
    figure(Name="Liss")
    p1 = plot(gammax, gammay);
    p1.LineWidth = 2;
    hold all
    % titel en
    % labels
    title('Lissajousfiguur en benaderende spline')
    xlabel('x') 
    ylabel('y')
    
    % respectievelijk de x- en y-waarden van de lissajousfiguur
    % interpoleren en evalueren over k punten
    % de eerste figuur met periodische voorwaarden en de tweede met natuurlijke
    [c_lisx, y_lisx] = cubicsplinesolve(t, tbefore, tafter, 1-i, [0; lisx; 0], k);
    [c_lisy, y_lisy] = cubicsplinesolve(t, tbefore, tafter, 1-i, [0; lisy; 0], k);
    
    % k geëvalueerde punten in de spline plotten
    p2 = plot(y_lisx, y_lisy, '--');
    p2.LineWidth = 2;
    % legende toevoegen
    legend('Lissajousfiguur', 'Spline', 'Location', 'northwest', 'FontSize', 12)
    hold off

end

% --- OPDRACHT 16 --- %

m = 2^13;
V = linspace(0, 2*pi, m);
[gammax, gammay] = lissajous(2, 3, pi/4, m);

tbefore = [-0.3*pi, -0.2*pi, -0.1*pi];
tafter = [2.1*pi, 2.2*pi, 2.3*pi];

interpolatiefout = zeros(1, 10);
r = 1:10;
amtPoints = zeros(1, length(r));
for i = r
    amtPoints(i) = 2^i;
end

for i = 1:length(r)
   t = linspace(0, 2*pi, amtPoints(i));
   [f_x, f_y] = lissajous(2, 3, pi/4, amtPoints(i));
   [c_lisx, y_lisx] = cubicsplinesolve(t, tbefore, tafter, 1, [0; f_x; 0], m);
   [c_lisy, y_lisy] = cubicsplinesolve(t, tbefore, tafter, 1, [0; f_y; 0], m);
   diffx = y_lisx - gammax.';
   diffy = y_lisy - gammay.';
   norms = sqrt(diffx.^2 + diffy.^2);
   interpolatiefout(i) = max(norms);
end

% uiteindelijk plotten we 2^r (het aantal knooppunten) en de interpolatiefout
figure(Name="Interpolatiefout bij de Lissajousfiguur")
semilogx(amtPoints, interpolatiefout)
title('Interpolatiefout bij de Lissajousfiguur')
xlabel('Aantal knooppunten') 
ylabel('Interpolatiefout')