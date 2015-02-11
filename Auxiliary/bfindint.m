function out = bfindint(A,s)

% Do a binary search for s through a sorted array A to find the i such that
% A(i) < s <= A(i+1). If A is empty, out is set to empty.

% Juan Kuntz, 06/02/2015

if isempty(A)
   out = [];
   return
end
if isscalar(A) % Do not actually need to compute exact length, just check whether it is less than 2 -- check whether we can do this better.
    if A < s
        out = 1;
        return
    else
        out = 0;
        return
    end
end

index = ceil(length(A)/2);
candidatel = A(index); candidater = A(index+1);

if candidatel < s && s <= candidater 
    out = index;
    return
elseif candidatel == s
    out = index - 1;
    return
elseif candidater + 1 == s
    out = index + 1;
    return
elseif candidater < s
    out = index + bfindint(A(index+1:end),s);
else
    out = bfindint(A(1:index-1),s);
end
    
end