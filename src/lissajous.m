function [gammax,gammay] = lissajous(a, b, phi, k)
% a, b en phi definiëren de parametrisatie van de lissajousfiguur
% gammax en gammay zijn vectoren van lengte k, die de x en y waarden geven voor theta % in linspace(0,2pi,k)

gammax = zeros(k, 1);
gammay = zeros(k, 1);

theta = 0;
increment = (2*pi)/(k-1);
for i = 1:k
    gammax(i) = sin(a*theta + phi);
    gammay(i) = sin(b*theta);
    theta = theta + increment;
end

end