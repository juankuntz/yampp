function out = iszero(p)

% Tests whether or not p is the zero polynomial, if it is iszero.m returns
% 1 otherwise 0.

% To do: overload so that it accepts matrix inputs.

% Juan Kuntz, 16/02/2015, last edited 16/02/2015.

out = 0;

if isempty(p.coef) || sum(p.coef(1,:)) == 0
   out = 1; 
end

end

