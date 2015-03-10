function obj = updatenvar(obj,n)

% Auxiliary function to set.var method that updates the total number of variables
% property of polynomial obj. A seperate function is required for this
% because set methods of some property are not allowed to set other
% properties of the object.

% Juan Kuntz, 05/02/2015, last edited 10/03/2015

obj.nvar = n;

end

