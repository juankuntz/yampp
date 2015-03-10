function out = decvar(obj,s,n)

% An auxiliary method called by /Auxiliary/vars.m to declare new
% independent variables.

% Juan Kuntz, 10/03/2015.

for i = 1:n
    temp.symb = s; temp.ncomp = n;
    obj(i,1).var = temp;
    obj(i,1).coef = [1;n+2-i]; % Declare coefficients in grlex order. However, this must be done after the variables are declared (otherwise the set.var function jumbles the coefficients).
    obj(i,1).deg = 1;
    obj(i,1).nvar = n;
    clear temp
end

out = obj;

end