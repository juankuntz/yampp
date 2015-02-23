function g = stimes(p,q)

% Computes the product of two scalars p and q.

% Juan Kuntz, 08/02/2015, last edited: 16/02/2015

% If either p or q are real numbers the product is easy.

if isdouble(p)
    g = q; 
    if isempty(q.coef) % If q is the zero polynomial we are done.
        return
    else % Otherwise just multiply each coefficient of q by p.
        g.coef(1,:) = p*q.coef(1,:);
    end
    return
end
if isdouble(q)
    g = stimes(q,p);
    return
end

% Make sure that p and q have the same variables and make template for g.

if strcmp(p.var.symb,q.var.symb)
    g = p; g.coef = []; 
else
    [pnotq,qnotp] = varcomp(p,q);
    if ~isempty(qnotp.symb)
        p.var = qnotp; 
    end
    if ~isempty(pnotq)
        q.var = pnotq;
    end
	g = p; g.coef = []; 
end

% Set degree of g and generate choosetables

g.deg = p.deg+q.deg; g.choose = [];
g.choose = ncktab(g.nvar+g.deg);

% Compute all individual products.

n = sum(p.var.ncomp);

if isempty(p.coef) || isempty(q.coef) % If either is the zero polynomial, return the zero polynomial.
    return
end

I = numel(p.coef(1,:)); J = numel(q.coef(1,:));

% Here we only need to call grlex I+J times and not IxJ times.

tempp = zeros(n,I); tempq = zeros(n,J);

for i = 1:I
    tempp(:,i) = grlext(n,p.coef(2,i),g.choose);
end

for j = 1:J
    tempq(:,j) = grlext(n,q.coef(2,j),g.choose);
end

% The product of each pair of monomials.

g.coef = zeros(2,J*I); l = 1;
for i = 1:I
    for j = 1:J
        g.coef(1,l) = p.coef(1,i)*q.coef(1,j);
        g.coef(2,l) = igrlext(tempp(:,i)+tempq(:,j),g.choose);
        l = l + 1;
    end
end
clear I J tempp tempq

% Sort out the entries in g.

[g.coef(2,:),I] = sort(g.coef(2,:));
temp = g.coef(1,:); g.coef(1,:) = temp(I);
clear I

% Combine entries that correspond to the same monomial (removing the
% entries of zero coefficients).

i = 2; l = numel(g.coef(1,:));
while i <= l
   if g.coef(2,i-1) == g.coef(2,i) 
       if g.coef(1,i-1) == -g.coef(1,i) 
           g.coef(:,i-1:i) = [];
           i = i - 1;
           l = l - 2;
       else
           g.coef(1,i-1) = g.coef(1,i-1) + g.coef(1,i);
           g.coef(:,i) = [];
           l = l - 1;
       end
   else
       i = i + 1;
   end
end

end

