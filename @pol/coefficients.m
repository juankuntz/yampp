function out = coefficients(varargin)

% Returns the coefficients of a polynomial. Takes in two arguments:

% 1. A pol p.
% 2. A string s containing either 'full' or 'sparse'. If s is not specified
% it is set by default to 'full'.

% If s is 'full' out is a vector containing all coefficients of p indexed by
% grlex order (we include those only of the monomials of degree less or 
% equal than p.deg). If s is 'sparse' we return two arrays, the first
% containing the coefficients, the second the respective ranks in grlex.

% Juan Kuntz, 17/02/2015, last edited 04/03/2015

p = varargin{1};

if nargin == 1
    s = 'full';
else
    s = varargin{2};
end

switch s
    case 'full'

        N = nchoosek(p.nvar+p.deg,p.deg);
        J = numel(p.coef(1,:));

        j = 1;
        out = zeros(N,1);
        
        for i = 1:N
            rank = p.coef(2,j);
            if i == rank
                out(i) = p.coef(1,j);
                if j == J
                    return
                else
                    j = j + 1;
                end
            else
                out(i) = 0;
            end
        end
    case 'sparse'
        out = p.coef;
end

end