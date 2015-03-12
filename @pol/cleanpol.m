function g = cleanpol(p)

% Function inspects p and looks for symbols that do not appear in any of
% p's monomials with zero coefficients, but that do appear in the list of
% variables of p. If it finds any, it removes them and outputs p rewritten
% in the smaller list of variables. If it doesn't, it just returns p.

% For matrix arguments, it looks for symbols that do not appear in any of
% the entries of the matrix.

% It also updates the degree in the case that the highest degree monomials
% have been deleted.

% Juan Kuntz, 11/03/2005, last edited 12/03/2005.

% First decide which symbols are actually present in the polynomial, the
% ith symbol is present if and only if v(i)=1.

[n,m] = size(p);

DEG = 0;
for i = 1:n
    for j = 1:m
        mon{i,j} = monvectors(p(i,j));
        DEG = max([DEG,sum(mon{i,j}(:,end))]);
    end
end

% DEG now contains the degree of the new matrix polynomial.

V = [];
for i = 1:n
    for j = 1:m
        if ~isempty(mon{i,j})
            V = [V,whichvar(p(i,j),mon{i,j})];
        end
    end
end

v = any(V,2); clear V

% If all symbols are present, there is no symbol cleaning to do, we just
% need to update the degrees.

if sum(v) == numel(p(1,1).var.symb)
    g = p;
    for i = 1:n
        for j = 1:m
            g(i,j).deg = DEG;
        end
    end
    return
end

% Otherwise declare g polynomial in new variables

u = find(v);

gtemp = pol;
temp.symb = p(1,1).var.symb(u);
temp.ncomp = p(1,1).var.ncomp(u);
gtemp.var = temp;

g(n,m) = gtemp;

clear temp gtemp

% Rewrite the polynomials in the reduced set of variables.

for i = 1:n
    for j = 1:m
        g(i,j) = rewrite(p(i,j),v,u,mon{i,j});
        g(i,j).deg = DEG;
    end
end
end


function v = whichvar(p,mon)

v = zeros(numel(p.var.symb),1);

s = 1;
for i = 1:numel(p.var.symb)
    v(i) = any(any(mon(s:s+p.var.ncomp(i)-1,:)));
    s = s + p.var.ncomp(i);
end

end


function g = rewrite(p,v,u,mon)

% Declare g polynomial in new variables

g = pol;

% If there are zero symbols in the end, we just need to need to decide
% whether the polynomial is a nonzero constant or just zero. Also set nvar
% to zero.

if isempty(u)
    if ~iszero(p)
       g.coef = p.coef; 
    end
    g.nvar = 0;
    return
end

% Otherwise we need to do a bit more work, starting with setting the
% correct variables for g. Also update nvar.

temp.symb = p(1,1).var.symb(u);
temp.ncomp = p(1,1).var.ncomp(u);
g.var = temp;
g.nvar = sum(p(1,1).var.ncomp(u));

clear temp

% Re-write the coefficients of p in terms of the grlex order over the
% reduced number of variables and store them in g.

L = sum(p.var.ncomp(u)); % Store these to avoid unnecessary computation.
TAB = ncktab(L+p.deg);  

for i = 1:numel(p.coef(1,:))
    newmon = zeros(L,1);
    l = 1; s = 1; k = 1;
    for j = 1:numel(p.var.symb)
       if v(j) == 1
           newmon(s:s+p.var.ncomp(u(l))-1) = mon(k:k+p.var.ncomp(j)-1,i);
           l = l + 1;
           s = s + p.var.ncomp(u(l-1));
       end
        k = k + p.var.ncomp(j);
    end
    g.coef(:,end+1) = [p.coef(1,i);igrlext(newmon,TAB)];
end

end