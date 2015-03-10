function out = isvar(x)

% Checks whether or not each entry of x is a variable in the sense that it
% is a first degree monomial with unity coefficient. In addition, we do not
% allow vectors with repeated entries count as vectors of variables. For
% example
%
% >> pol('x',2); isvar(x)
%
% returns 1, while
%
% >> pol('x',2); isvar([x(1),x(1)])
%
% returns 0. Currently, only supports vectors of polynomials.

% Juan Kuntz, 02/02/2015, last edited 10/03/2015

if min(size(x)) ~= 1
   error('x cannot be a matrix');
end

% If x has any repeated entries, then it is not a vector of variables (our
% definition).

for i = 1:numel(x)
    for j = i+1:numel(x)
        if x(i) == x(j)
            out = 0;
            return
        end
    end
end

% If any entry of x has more than one nonzero monomial, if it is not of
% degree one, or if it does not have unity coefficient, then x is not a
% vector of variables.

for i = numel(x)
    if ~isscalar(x(i).coef(1,:)) || x(i).coef(1,1) ~= 1 || nnz(igrlext(x(i).coef(2,1),ncktab(x(i).nvar+x(i).deg))) ~= 1 
       out = 0; 
       return
    end
end

% If we've gotten this far, then x is a vector of variables!

out = 1;
end