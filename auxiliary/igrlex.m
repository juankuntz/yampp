function out = igrlex(a)

% Returns the grlex rank of the multiindex a. Essentially just reverts what
% the function grlex.m does. 

% Juan Kuntz, 06/02/2015

n = numel(a); d = sum(a);

if d == 0
    out = 1;
    return
end

out = binomial(n+d-1,d-1);

for i = 1:n
   out = out + binomial(n-i+d,d) - binomial(n-i+d-a(i),d-a(i));
   d = d - a(i);
end

out = out + 1; % Add 1 because we rank the zero monomial as 1 in MATLAB 
% instead of 0.
end