function out = eq(p,q)

% Returns 1 if p and q are the same polynomials. If they are not scalar, eq
% returns 1 if p(i,j) and q(i,j) are the same polynomials for all i,j.
% Otherwise returns 0.

% Juan Kuntz, 10/03/2015.

out = 1;

[n,m] = size(p);

% If p and q do not have the same dimension, then they are not equal.

if ~isequal([n,m],size(q)) 
    out = 0;
    return
end

% If p and q are not scalars, test equality elementwise.

if max(size(p))>1
    for i = 1:n
        for j = 1:m
            if eq(p(i,j),q(i,j)) == 0
                out = 0;
                return
            end
        end
    end
end

% If p and q are not in the same variables, then they are not equal.

if ~isequal(p.var.symb,q.var.symb) || ~isequal(p.var.ncomp,q.var.ncomp) 
    out = 0;
    return
end

ncoef = numel(p.coef(1,:));

% If p and q do not have the same number of non-zero coefficients, then
% they are not equal.

if ncoef ~= numel(q.coef(1,:));
    out = 0;
    return
end

% If p and q do not have the non-zero coefficients, and corresponding
% monomials, then they are not the same.

for i = 1:ncoef
    if p.coef(1,i) ~= q.coef(1,i) || p.coef(2,i) ~= q.coef(2,i)
       out = 0;
       return
    end
end

% If we've reached this far, then p and q are the same polynomials!

end

