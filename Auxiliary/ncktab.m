function out = ncktab(n)

% Creates table out such that out(i,j) = i-1 choose j-1 for all i
% <= n+1 and j <= i.

% Juan Kuntz 08/02/2015

out = zeros(n);

for i = 1:n+1
    for j = 1:i
        out(i,j) = nchoosek(i-1,j-1);
    end
end
end