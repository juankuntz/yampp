function out = ismon(p)

% Check whether or not p is a monomial.

% Juan Kuntz, 09/02/2015

if numel(p.coef(1,:)) > 1 || p.coef(1,1) ~= 1
    out = 0;
else 
    out = 1;
end


end

