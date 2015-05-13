function g = plus(p,q)

% Computes the componentwise sum between two polynomials p and q of the
% same dimensions.

% Juan Kuntz, 08/02/2015, last edited 12/05/2015

[n,m] = size(p); 

if ~isequal([n,m],size(q))
    error('Both polynomials to be summed must have the same dimension as vectors ');
end

flg = 0;
for i = 1:n
    for j = 1:m
        [g(i,j),temp] = splus(p(i,j),q(i,j));
        flg = max(temp,flg);
    end
end

% Some monomials have cancelled out, hence we may have to clean.

if flg
    g = cleanpol(g); % Removes uneeded symbols and ensures that degree is what it should be.
end
end

function [g,cleanflag] = splus(p,q)

% Computes the sum of two scalars. The cleanflag warns us that we should
% run clean pol later since at least some monomial has been canceled out
% and thus we may have removed some variables/lowered the degree.

cleanflag = 0;

% If p or q is a double.

if isdouble(p) 
    g = q;
    if p~=0
        if g.coef(2,1) == 1 % q has a nonzero zero-monomial term.
            g.coef(1,1) = g.coef(1,1) + p;
        else 
            g.coef = [[p;1],g.coef];
        end
    end
    return
end
if isdouble(q) 
    g = splus(q,p);
    return
end

% Make sure that p and q are written in the same variables and set g to be
% p.

if ~strcmp(p.var.symb,q.var.symb)
    [pnotq,qnotp] = varcomp(p,q);
    if ~isempty(qnotp.symb)
        p.var = qnotp; 
    end
    if ~isempty(pnotq.symb)
        q.var = pnotq;
    end
end

g = p;

% If q or p are the zero polynomial this is easy.

if isempty(q.coef)
    return
elseif isempty(p.coef)
    clear g
    g = q;
    return
end

% Set degree of g 

g.deg = max(p.deg,q.deg);

% Add q into g.

i = 1; j = 1; nc = numel(q.coef(1,:)); l = numel(g.coef(1,:));
while i <= l
    if g.coef(2,i) == q.coef(2,j)
        if g.coef(1,i) + q.coef(1,j) == 0 % We do not store zero coeficients
            cleanflag = 1;  % We may have to clean.
            g.coef(:,i) = [];
            l = l - 1;
        else
            g.coef(1,i) = g.coef(1,i) + q.coef(1,j); % New coeficient is the sum of the previous.
            i = i + 1;
        end
        j = j + 1;
    elseif g.coef(2,i) > q.coef(2,j)
        g.coef = [g.coef(:,1:i-1),q.coef(:,j),g.coef(:,i:end)];
        j = j + 1;
        i = i + 1;
        l = l + 1;
    else
        i = i + 1;
    end
    if j > nc
        break
    end
end

g.coef = [g.coef,q.coef(:,j:end)];

end