function out = grlext(varargin)

% Takes in two arguments, first n, then rank. Returns the multiindex
% rankend rank in the graded lexicographic order of monomials in n
% variables.

% Juan Kuntz, 05/02/2015

n = varargin{1}; rank = varargin{2};

% If there is only one variable, this is trivial. Or if the rank is 1, then
% this is also trivial.

if n == 1
   out = rank - 1; 
   return
elseif rank == 1 && nargin == 3 
    out = zeros(n,1);
    return
end

% Check whether we are supplied with the matrix full of binomial
% coefficients of the appropiate size (improves speed). If we are, then we
% don't have to bother calling repeatedly the function binomial.

% Remember the entries of tab are shifted by one, so tab(i,j) is i-1 choose
% j-1.

if nargin == 3 
    tab = varargin{3};
else
    tab = varargin{4};
end

% Compute the degree of the monomial: 
%
% - The first binomial(n+0,0) ranks correspond to the monomials of degree
% 0.
% - The next binomial(n+1,1) ranks correspond to the monomials of degree 1.
% - etc...

if nargin == 2 || numel(varargin{3}(1,:))>1 % We only need to compute the degree once.
    d = 1;
    while ~(tab(n+d,d)<rank && rank <= tab(n+d+1,d+1))
        d = d + 1;
    end
    rank = rank - tab(n+d,d);
    out = grlext(n,rank,d,tab); % Once we have figured out the degree of the 
    % monomials we can start figuring out exactly which monomial of that
    % degree the degree corresponds to. Note we remove the bit of the rank
    % that corresponds to the monomials of lower degree.
    return
end

dp = varargin{3};

for d = dp:-1:0 % Multiindexes ranked binomial(n-1+dp,dp)-binomial(n-1+d,d)
    % are those that have dp-d in this entry and the remaining entries sum up to d.
    if d == 0 || (tab(n+dp,dp+1)-tab(n+d,d+1) < rank && rank <= tab(n+dp,dp+1)-tab(n+(d-1),d))
        break
    end
end

if n == 2 % When there's only two entries to figure out, there's no more ambiguity.
    out = [dp-d;d];
    return
end

rank = rank -  (tab(n+dp,dp+1)-tab(n+d,d+1)); 
out = [dp-d;grlext(n-1,rank,d,tab)]; % Figure out what goes in the remaining entries of the multiindex.
end