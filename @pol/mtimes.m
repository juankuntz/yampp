function g = mtimes(p,q)

% Compute the matrix product of p and q (or the appropiate scalar product
% at least one of them is scalar. 

% Juan Kuntz, 08/02/2015, last edited 12/03/2015

% Case 1: Inner product between sdpvar (or seq) object and a pol object.

if isa(q,'sdpvar') || isa(q,'seq')
    if (min(size(q)) ~= 1 && isa(q,'sdpvar')) || ~isscalar(p) || (~isscalar(q) && isa(q,'seq'))
        disp('Error: Inner product between sdpvar and pol classes only implemented for scalars.');
        return
    end
    
    % Compute inner product.
    
    g = 0;
    
    if isempty(p.coef) % If p is the zero polynomial, return zero.
        return
    end
    
    np = numel(p.coef(1,:));
    j = 1;
    for i = 1:numel(q)
        if j <= np && p.coef(2,j) == i
            g = g + p.coef(1,j)*q(i);
            j = j + 1;
        end
    end
    return
end

% Case 2: Normal scalar or matrix multiplcation between pols and/or a pol
% and a real number.

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
  
% Now compute product.
    
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
    error(['The number of columns of ',p,' must be the same as that of rows of ',q, ' in order for their matrix product to make sense, otherwise one of them must be a scalar to compute the scalar product']);
end

if cflg
    g = cleanpol(g);
end
end

