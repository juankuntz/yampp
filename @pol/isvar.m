function out = isvar(x)

% Checks whether or not each entry of x is a variable in the sense that it
% is a first degree monomial with unity coefficient. Currently only works
% for vectors of polynomials.

% Juan Kuntz, 02/02/2015, last edited 10/03/2015

if min(size(x)) ~= 1
   disp('Error: x cannot be a matrix');
   return
end

out = 1;

for i = numel(x)
    if ~isscalar(x(i).var.symb) || ~isscalar(x(i).coef(1,:)) || x(i).coef(1,1) ~= 1 || nnz(igrlext(x(i).coef(2,1),ncktab(x(i).nvar+x(i).deg))) ~= 1 
       out = 0; 
       return
    end
end

end