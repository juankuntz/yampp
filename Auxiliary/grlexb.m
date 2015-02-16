function out = grlexb(n,d)

% Returns the multiindexes of degree d in n variables order in the grlex
% order.

% Juan Kuntz, 06/02/2015

out = [];

% Easy base cases

if n == 1
    out = d;
    return
elseif n == 2
    out = [0:d;d:-1:0];
    return
end

% Otherwise recurse

for i = 1:d+1
    temp = [(i-1)*ones(1,nchoosek(n-1+d-(i-1)-1,d-(i-1))); grlexb(n-1,d-(i-1))];
    out = [out,temp];
    clear temp
end
end