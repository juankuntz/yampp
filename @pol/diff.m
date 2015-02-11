function dp = diff(p,x)

% Computes the jacobian [dp(1)/dx(1),dp(1)/dx(2),...;
% dp(2)/dx(1),dp(2)/dx(2),...; ... ].

% Juan Kuntz, 09/02/2015.

% Check that x is a vector of variables.

if ~isvar(x) 
    disp(['Error: You can only differentiate with respect to a scalar variable (or a vector of scalar varaibles); ',x,' is neither']);
    return
end

% Check that p is not a matrix.

if min(size(p)) > 1
    disp('Error: Differentiation of a matrix of polynomials is not supported in this version.')
    return
end

% Compute derivatives one by one.

for i = 1:numel(p)
    for j = 1:numel(x)
        dp(i,j) = sdiff(p(i),x(j));
    end
end


end

function dp = sdiff(p,x)

% There is a smart way of doing this that does not require calling
% grlex/iglrex (there's clearly a pattern in how the ranks change when we
% differentiate, just try it out.

[pnotx,xnotp] = varcomp(p,x);

if ~isempty(xnotp.symb)
    dp = pol;
    return
elseif ~isempty(pnotx.symb)
    x.var = pnotx;
end

monox = grlext(x.nvar,x.coef(2,1),x.choose); 
dp = p; dp.coef = [];

for i = 1:numel(p.coef(1,:))
    temp = grlext(p.nvar,p.coef(2,i),p.choose);
    if temp(logical(monox)) ~= 0 % That is, if variable x appears in the monomial.
        dp.coef = [dp.coef,[temp(logical(monox))*p.coef(1,i);igrlext(temp-monox,x.choose)]];
    end
end
clear temp

if ~isempty(dp.coef) % Make sure that monomials are in the correct order. 
    [temp,I] = sort(dp.coef(2,:));
    dp.coef = [dp.coef(1,I);temp];
end
end