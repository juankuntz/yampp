function g = vertcat(varargin)

% Overloads vertcat so that the variables of the input matrices of
% polynomials are updated, for example g = [p,q] returns the concatenation
% of x and y, but the variables of all entries of p include both those p
% and those of q. This is essentially just a copy of the code for horzcat,
% with a couple of tranposes.

% Juan Kuntz, 11/03/2015.

k = 0;
for i = 1:nargin
    if ~isempty(varargin{i})
        k = k + 1;
        temp{k} = varargin{i}';
    end
end

for i = 1:k
    [n(i),m(i)] = size(temp{i});
end

for i = 2:k
    if n(i-1) ~= n(i)
        error('To concatenate two matrices of polynomials horizontally, each matrix must have the same number of rows');
    end
end

% Initialise g and adding the variables of temp{1} into it.

g(n(1),sum(m)) = pol;

% Populate g

for i = 1:k
    
    p = temp{i};
    
    % Make sure that g and p have the same variables.
    
    if ~strcmp(g(1,1).var.symb,p(1,1).var.symb)
        [gnotp,pnotg] = varcomp(g(1,1),p(1,1));  
        for j = 1:n(1)
            for k = 1:sum(m)
                if ~isempty(pnotg.symb)
                    g(j,k).var = pnotg; 
                end
            end
            for k = 1:m(i)
                if ~isempty(gnotp.symb)
                    p(j,k).var = gnotp;
                end
            end
        end
    end

    g(:,sum(m(1:i-1))+1:sum(m(1:i))) = p; % Add in the next matrix.
    
    clear p gnotp pnotg
end

TEMP = g; clear g; g = TEMP';

end