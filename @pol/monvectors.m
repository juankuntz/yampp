function out = monvectors(p)

% Returns an array containing the multiindexes of the monomials with
% nonzero coefficients in p, as column vectors, and in the order in which
% the ranks of the monomials are stored in p.coef.

% Juan Kuntz, 11/03/2015.

if isempty(p.coef)  % Zero polynomial
    out = zeros(p(1,1).nvar);
    return
end
            
TAB = ncktab(p.nvar+p.deg);
l = numel(p.coef(1,:));

out = zeros(p.nvar,l);

for i = 1:l
    out(:,i) = grlext(p.nvar,p.coef(2,i),TAB);
end

end