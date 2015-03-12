function g = times(p,q)

% Compute the elementwise product of p and q (or the appropiate scalar product
% at least one of them is scalar.

% Juan Kuntz, 09/02/2015, last edited 12/03/2015

[np,mp] = size(p); [nq,mq] = size(q);

% If we are multiplying by zeros we have to clean later in case we've
% removed some symbol or lowered the degree.

cflg = 0;
if isa(p,'double')
    cflg = max(cflg,any(any(p==0)));
    cflg = max(cflg,any(any(iszero(q))));
elseif isa(q,'double')
    cflg = max(cflg,any(any(q==0)));
    cflg = max(cflg,any(any(iszero(p))));
else
    cflg = max(cflg,any(any(iszero(q))));
    cflg = max(cflg,any(any(iszero(p))));
end

if np == 1 && mp == 1 % Scalar product
    for i = 1:nq
        for j = 1:mq
            g(i,j) = stimes(p,q(i,j));
        end
    end
elseif nq == 1 && mq == 1 % Scalar product
    for i = 1:np
        for j = 1:mp
            g(i,j) = stimes(q,p(i,j));
        end
    end
elseif isequal([np,mp],[nq,mq]) % Elementwise product
    for i = 1:np
        for j = 1:mp
            g(i,j) = stimes(p(i,j),q(i,j));
        end
    end
else
    error(['Error: ',p,' must have the same dimensions as ',q, ' in order for their elementwise product to make sense, otherwise one of them must be a scalar to compute the scalar product.'])
end

if cflg
    g = cleanpol(g);
end

end

