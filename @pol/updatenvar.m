function [obj,oldchoose] = updatenvar(obj,n)

% Auxiliary function to set.var that updates the total number of variables
% property of polynomial obj.

obj.nvar = n;
oldchoose = obj.choose;
obj.choose = []; obj.choose = ncktab(n+obj.deg); 
end

