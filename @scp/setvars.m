function setvars(cp)

% Decide what's the dimension of the underlying space. Currently, only
% polynomials in a single vector of variables are supported.

% Juan Kuntz, 13/03/2015.

p = cp.obj{2}(1);

cp.nvar = prop(p,'n');

end

