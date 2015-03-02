function v = monlist(x,d)

% Returns the vector of monomials in variables x of degree d or less.

% Juan Kuntz, 02/03/2015

n = numel(x);

tab = ncktab(n+d);

v = pol(1);
for i = 2:nchoosek(n+d,d)
    mon = grlext(n,i,tab);
    v(i) = x^mon;
    clear mon
end

end

