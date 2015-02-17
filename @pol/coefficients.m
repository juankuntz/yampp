function out = coefficients(varargin)

% Returns the coefficients of a polynomial. Takes in two arguments:

% 1. A pol p.
% 2. A string s containing either 'full' or 'sparse'. If s is not specified
% it is set by default to 'full'.

% If s is 'full' out is a vector containing all coefficients of p indexed by
% grlex order (we include those only of the monomials of degree less or 
% equal than p.deg). If s is 'sparse' we return two arrays, the first
% containing the coefficients, the second the respective ranks in grlex.

% Juan Kuntz, 17/02/2015, last edited 17/02/2015

p = varargin{1};

if nargin == 1
    s = 'full';
else 
    s = varargin{2};
end

switch s
    case 'full'
        j = 1;
        out = zeros(nchoosek(p.nvar+p.deg,p.deg),1);
        for i = 1:nchoosek(p.nvar+p.deg,p.deg)
            if i == p.coef(2,j)
                out(i) = p.coef(1,j);
                if j == numel(p.coef(1,:))
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