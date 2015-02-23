function out = ncktab(n)

% Creates table out such that out(i,j) = i-1 choose j-1 for all i
% <= n+1 and j <= i. This table is stored for future use in the global
% variable NCKTAB. If NCKTAB already exists and is at least as large as
% (n+1)^2, then we just return out = NCKTAB. We do this to speed up
% computations by avoiding unnecessary calls of nchoosek.

% Juan Kuntz 08/02/2015, last edited 23/02/2015

global NCKTAB

if isempty(NCKTAB) || NCKTAB.n < n
    NCKTAB = [];
    NCKTAB.tab = zeros(n);
    NCKTAB.n = n;
    for i = 1:n+1
        for j = 1:i
            NCKTAB.tab(i,j) = nchoosek(i-1,j-1);
        end
    end
end

out = NCKTAB.tab;

end


