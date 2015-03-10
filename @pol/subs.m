function g = subs(p,q)

% Returns g = p(q(1),q(2),...,q(end)).

% Juan Kuntz, 09/02/2015

if numel(p.var.symb)>1
    disp(['Error: ',p,' has variables in more than one symbol, the current version of subs only supports polynomials "that are being subbed into" which have a single symbol']);
    return
end

if min(size(q))>1 || numel(q) ~= p.nvar
    disp(['Error: ',q,' must have the same number of components as ',p,' has of variables']);
    return
end

% Generate template for g

g = q(1); % g has the same variables as q.
g.coef = []; % g is the zero polynomial for now.

TAB = ncktab(p.nvar+p.deg);
for i = 1:numel(p.coef(1,:))
    monop = grlext(p.nvar,p.coef(2,i),TAB);
    g = g + p.coef(1,i)*q^monop;
end

end