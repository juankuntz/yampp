function out = bfind(A,s)

% Do a binary search for s through a sorted array A. Returns empty argument
% if s is not in array A.

% Juan Kuntz, 06/02/2015

if isempty(A)
    out = [];
    return
end

index = ceil(numel(A)/2);
candidate = A(index);

if candidate == s
    out = index;
    return
elseif candidate < s
    out = index + bfind(A(index+1:end),s);
else
    out = bfind(A(1:index-1),s);
end
    
end