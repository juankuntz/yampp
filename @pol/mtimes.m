function g = mtimes(p,q)

% Compute the matrix product of p and q (or the appropiate scalar product
% at least one of them is scalar.

% Juan Kuntz 08/02/2015

[np,mp] = size(p); [nq,mq] = size(q);

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
elseif mp == nq % Matrix product
    for i = 1:np
        for j = 1:mq
            g(i,j) = stimes(p(i,1),q(1,j));
            for k = 2:mp
                g(i,j) = g(i,j) + stimes(p(i,k),q(k,j));
            end
        end
    end
else
    disp(['Error: the number of columns of ',p,' must be the same as that of rows of ',q, ' in order for their matrix product to make sense, otherwise one of them must be a scalar to compute the scalar product'])
end
end

