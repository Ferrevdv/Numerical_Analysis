function V=Vandermonde(x)
% x is een verticale vector met n verschillende elementen
% V is de bijhorende nxn Vandermonde-matrix

n = length(x);
V = ones(n,n);
for i = 1:n
    for j = 2:n     % de eerste kolom is steeds gevuld met enen.
        V(i, j) = x(i)^(j-1);
    end
end