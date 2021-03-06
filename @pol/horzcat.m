function g = horzcat(varargin)

% Overloads horzcat so that the variables of the input matrices of
% polynomials are updated, for example g = [p,q] returns the concatenation
% of x and y, but the variables of all entries of p include both those p
% and those of q.

% Juan Kuntz, 10/03/2015, last edited 12/03/2015.

% Remove empty polynomials and compute new degree

k = 0; DEG = 0;
for i = 1:nargin
    if ~isempty(varargin{i})
        k = k + 1;
        if isa(varargin{i},'double')
            varargin{i} = pol(varargin{i});
        end
        temp{k} = varargin{i};
        DEG = max(DEG,temp{k}(1,1).deg);
    end
end

for i = 1:k
    [n(i),m(i)] = size(temp{i});
end

for i = 2:k
    if n(i-1) ~= n(i)
        error('To concatenate two matrices of polynomials horizontally, each matrix must have the same number of rows.');
    end
end

% Initialise g.

g(n(1),sum(m)) = pol;

% Populate g and set the correct variables.

for i = 1:k
    
    p = temp{i};
    
    % Make sure that g and p have the same variables.
    
    if ~strcmp(g(1,1).var.symb,p(1,1).var.symb)
        [gnotp,pnotg] = varcomp(g(1,1),p(1,1));  
        if ~isempty(pnotg.symb)
            for j = 1:n(1)
                for k = 1:sum(m)
                    g(j,k).var = pnotg; 
                end
            end
        end
        
        if ~isempty(gnotp.symb)
            for j = 1:n(1)
                for k = 1:m(i)
                    p(j,k).var = gnotp;
                end
            end
        end
    end

    g(:,sum(m(1:i-1))+1:sum(m(1:i))) = p; % Add in the next matrix.
    
    clear p gnotp pnotg
end

% Update the degrees.

for i = 1:n(1)
    for j = 1:sum(m)
        g(i,j).deg = DEG;
    end
end
    
end