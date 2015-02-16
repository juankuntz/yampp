function out = sum(varargin)

% An overloaded version of the standard sum.m that takes in matrices of pol 
% objects.

% Juan Kuntz, 16/02/2015

A = varargin{1};
[n,m] = size(A);

% If A is a vector just add up the entries.

if min(n,m) == 1
    out = A(1);
    for i = 2:max(n,m)
        out = out + A(i);
    end
    return
end

% Otherwise, sum up the elements in each column.

if nargin == 1 || varargin{2} == 1
    for j = 1:m
        out(1,j) = A(1,j);
        for i = 2:n
            out(1,j) = out(1,j) + A(i,j);
        end
    end
    return
end


end