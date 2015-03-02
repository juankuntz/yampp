function v = monlist(x,d)

% Returns the vector of monomials in variables x of degree d or less.

% Juan Kuntz, 02/03/2015

n = numel(x);

tab = nchoosek(n+d);

for i = 1:nchoosek(n+d,d)
    mon = igrlext(i,tab);
    v(i) = x^mon;
    clear mon
end

end

