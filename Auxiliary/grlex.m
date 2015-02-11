function out = grlex(varargin)

% Takes in two arguments, first n, then rank. Returns the multiindex
% rankend rank in the graded lexicographic order of monomials in n
% variables.

% Juan Kuntz 05/02/2015

n = varargin{1}; rank = varargin{2};

% If there is only one variable, this is trivial.

if n == 1
   out = rank - 1; 
   return
elseif rank == 1
    out = zeros(n,1);
    return
end

% Compute the degree of the monomial: 
%
% - The first nchoosek(n+0,0) ranks correspond to the monomials of degree
% 0.
% - The next nchoosek(n+1,1) ranks correspond to the monomials of degree 1.
% - etc...

if length(varargin) == 2 % We only need to compute the degree once.
    d = 1;
    while ~(nchoosek(n+d-1,d-1)<rank && rank <= nchoosek(n+d,d))
        d = d + 1;
    end
    rank = rank - nchoosek(n+d-1,d-1);
    out = grlex(n,rank,d); % Once we have figured out the degree of the 
    % monomials we can start figuring out exactly which monomial of that
    % degree the degree corresponds to. Note we remove the bit of the rank
    % that corresponds to the monomials of lower degree.
    return
end

dp = varargin{3};

for d = dp:-1:0 % Multiindexes ranked nchoosek(n-1+dp,dp)-nchoosek(n-1+d,d)
    % are those that have dp-d in this entry and the remaining entries sum up to d.
    if d == 0 || (nchoosek(n-1+dp,dp)-nchoosek(n-1+d,d) < rank && rank <= nchoosek(n-1+dp,dp)-nchoosek(n-1+(d-1),d-1))
        break
    end
end

if n == 2 % When there's only two entries to figure out, there's no more ambiguity.
    out = [dp-d;d];
    return
end

rank = rank -  (nchoosek(n-1+dp,dp)-nchoosek(n-1+d,d)); 
out = [dp-d;grlex(n-1,rank,d)]; % Figure out what goes in the remaining entries of the multiindex.
end