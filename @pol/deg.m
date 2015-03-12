function d = deg(p)

% Returns the degree of the polynomial p, entrywise.

% Juan Kuntz, 12/03/2015.

[n,m] = size(p);

d = zeros(n,m);

TAB = ncktab(p(1,1).nvar+p(1,1).deg);

for i = 1:n
    for j = 1:m
        if ~isempty(p(i,j).coef)
            d(i,j) = sum(grlext(p(1,1).nvar,p(i,j).coef(2,end),TAB));
        end
    end
end

end