function out = iszero(p)

% Tests whether or not p is the zero polynomial, if it is iszero.m returns
% 1 otherwise 0.

% Juan Kuntz, 16/02/2015, last edited 13/03/2015.

[n,m] = size(p);

out = zeros([n,m]);

for i = 1:n
    for j = 1:m
        if isempty(p(i,j).coef) 
           out(i,j) = 1; 
        end
    end
end

end

