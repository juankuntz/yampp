function g = mpower(p,n)

% If n is a scalar, then g is the nth matrix power of p. If n is a
% multiindex then g = (p(1)^n(1))*(p(2)^n(2))*...

% Juan Kuntz, 08/02/2015, last edited 16/04/2015

if isscalar(n)
    g = pol(1);
    if n ~= 0
        for i = 1:n
            g = g*p;
        end
    end
    return
end

if min(size(p)) > 1 || numel(n) ~= numel(p)
    disp(['Error: ',n,' must be a scalar or a multiindex of the same length as the vector ',p]);
    return
end

g = 1;
for i = 1:numel(n)
    g = g*p(i)^n(i);
end

end


