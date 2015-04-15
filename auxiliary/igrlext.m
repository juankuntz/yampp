function out = igrlext(a,tab)

% Returns the grlex rank of the multiindex a. Essentially just reverts what
% the function grlex.m does. 

% Juan Kuntz, 06/02/2015

n = numel(a); d = sum(a);

if d == 0
    out = 1;
    return
end

out = tab(n+d,d);

for i = 1:n
   out = out + tab(n-i+d+1,d+1) - tab(n-i+d-a(i)+1,d-a(i)+1);
   d = d - a(i);
end

out = out + 1; % Add 1 because we rank the zero monomial as 1 in MATLAB 
% instead of 0.
end