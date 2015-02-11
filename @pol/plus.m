function g = plus(p,q)

% Computes the componentwise sum between two polynomials p and q of the
% same dimensions.

% Juan Kuntz, 08/02/2015

[n,m] = size(p); 

if ~isequal([n,m],size(q))
    disp(['Error: ',p,' and ',q, 'must have the same dimension as vectors in order to compute their sum']);
    return
end

for i = 1:n
    for j = 1:m
        g(i,j) = splus(p(i),q(i));
    end
end
end

function g = splus(p,q)

% Computes the sum of two scalars.

% If p or q is a real number this is easy

if isdouble(p) || isempty(p.coef) || numel(p.coef(1,:)) == 1 &&  p.coef(2,1) == 1
    g = q;
    if g.coef(2,1) == 1 % q has a nonzero zero-monomial term.
        if isdouble(p)
            g.coef(1,1) = g.coef(1,1) + p;
        elseif ~isempty(p.coef)
            g.coef(1,1) = g.coef(1,1) + p.coef(1,1);
        end
    else 
        if isdouble(p)
            g.coef = [[p;1],g.coef];
        elseif isempty(p.coef)
            g.coef = [p.coef,g.coef];
        end
    end
    return
end
if isdouble(q) || isempty(q.coef) || numel(q.coef(1,:)) == 1 &&  q.coef(2,1) == 1
    g = splus(q,p);
    return
end

% Make sure that p and q are written in the same variables and set g to be
% p.

if strcmp(p.var.symb,q.var.symb)
    g = p; 
else
    [pnotq,qnotp] = varcomp(p,q);
    if ~isempty(qnotp.symb)
        p.var = qnotp; 
    end
    if ~isempty(pnotq)
        q.var = pnotq;
    end
	g = p; 
end

% Set degree of g and generate choosetables

g.deg = max(p.deg,q.deg);
g.choose = ncktab(g.nvar+g.deg);

% Add q into g.

i = 1; j = 1; nc = numel(q.coef(1,:)); l = numel(g.coef(1,:));
while i <= l
    if g.coef(2,i) == q.coef(2,j)
        if g.coef(1,i) + q.coef(1,j) == 0 % We do not store zero coeficients
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
    else
        i = i + 1;
    end
    if j > nc
        break
    end
end

g.coef = [g.coef,q.coef(:,j:end)];

end