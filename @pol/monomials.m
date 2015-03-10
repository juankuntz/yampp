function v = monomials(x,d)

% Returns the vector of monomials, ordered in the grlex order, of degrees
% d(1), d(2), ...,  d(numel(d)) in variables x. It ignores repeated entries
% in d (because grlexb, the function used to generate the multiindexes),
% does so as well.

% Juan Kuntz, 10/03/2015.

if ~isvar(x)
    error('To generate the vector of monomials, the vector specifying the vector of monomials must be a vector of variables, that is a vector of monomials of degree 1, with unity coefficients, and no repeated entries');
end

inds = grlexb(numel(x),d); % Generate array of multiindexes.

v(numel(inds(1,:))) = pol; % Initialise the array of pols.

% Generate the monomials.

for i = 1:numel(inds(1,:))
    v(i) = x^inds(:,i);
end

end