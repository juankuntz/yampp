function dp = diff(p,x)

% Computes the jacobian [dp(1)/dx(1),dp(1)/dx(2),...;
% dp(2)/dx(1),dp(2)/dx(2),...; ... ].

% Juan Kuntz, 09/02/2015, last edited 11/03/2015.

% Check that x is a vector of variables.

if ~isvar(x) 
    error(['You can only differentiate with respect to a scalar variable (or a vector of scalar varaibles); ',x,' is neither']);
end

% Check that p is not a matrix.

if min(size(p)) > 1
    error('Differentiation of a matrix of polynomials is not supported in this version.')
end

% Compute derivatives one by one.

for i = 1:numel(p)
    for j = 1:numel(x)
        dp(i,j) = sdiff(p(i),x(j));
    end
end

if numel(dp(1,1).var.symb)> 1    % Only need to clean here if there is more than one symbol, otherwise we can use the cheaper code below (just need to check whether we've dropped degree/and or become a constant).
    dp = cleanpol(dp);
else
    DEG = 0;
    
    % Compute degree
    
    for i = 1:numel(p)
        for j = 1:numel(x)
            DEG = max(DEG,dp(i,j).deg);
        end
    end
    
    % Update degree
    
    for i = 1:numel(p)
        for j = 1:numel(x)
            dp(i,j).deg = DEG;
        end
    end
    
    % If dp is a constant, get rid of it's variables.
    
    if DEG == 0         
       temp = dp; clear dp; 
       dp = pol(zeros(size(temp))); 
       for i = 1:numel(p)
            for j = 1:numel(x)
                dp(i,j).coef = temp(i,j).coef;
            end
        end
    end
end

end

function dp = sdiff(p,x)

% There is a smart way of doing this that does not require calling
% grlex/iglrex (there's clearly a pattern in how the ranks change when we
% differentiate, just try it out.

if isdouble(p) % If p is a double return the zero polynomial.
    dp = x;
    dp.coef = [];
    return
end

[pnotx,xnotp] = varcomp(p,x);

if ~isempty(xnotp.symb) 
    dp = pol;
    return
elseif ~isempty(pnotx.symb)
    x.var = pnotx;
end

monox = grlext(x.nvar,x.coef(2,1),ncktab(x.nvar+x.deg)); 
dp = p; dp.coef = [];

if isempty(p.coef) % If p is the zero polynomial, return the zero polynomial.
    return
end

TAB = ncktab(p.nvar+p.deg);
for i = 1:numel(p.coef(1,:))
    temp = grlext(p.nvar,p.coef(2,i),TAB);
    if temp(logical(monox)) ~= 0 % That is, if variable x appears in the monomial.
        dp.coef = [dp.coef,[temp(logical(monox))*p.coef(1,i);igrlext(temp-monox,TAB)]];
    end
end
clear temp

if ~isempty(dp.coef) % Make sure that monomials are in the correct order and update the degree.
    [temp,I] = sort(dp.coef(2,:));
    dp.coef = [dp.coef(1,I);temp];
    dp.deg = sum(grlext(dp.nvar,dp.coef(2,end),ncktab(dp.deg+dp.nvar))); % Update degree.
else
    dp.deg = 0; % Update degree.
end

end