function out = prop(varargin)

% Used to check out the properties of pol p. We need since we've overloaded subsref.

% Juan Kuntz, 12/03/2015

p = varargin{1}; [n,m] = size(p);

if nargin == 1
    s = 'a';
else
    s = varargin{2};
end

switch s
    case 'c'
        if n == 1 && m == 1
            out = p.coef;
            return
        end
        
        for i = 1:n
            for j = 1:m
                out{i,j} = p(i,j).coef;
            end
        end
        
    case 'v'
        
        if n == 1 && m == 1
            out = p.var;
            return
        end
        
        for i = 1:n
            for j = 1:m
                out{i,j} = p(i,j).var;
            end
        end
        
    case 'n'
        
        if n == 1 && m == 1
            out = p.nvar;
            return
        end
        
        
        for i = 1:n
            for j = 1:m
                out{i,j} = p(i,j).nvar;
            end
        end
        
    case 'd'
        
        if n == 1 && m == 1
            out = p.deg;
            return
        end
        
        
        for i = 1:n
            for j = 1:m
                out{i,j} = p(i,j).deg;
            end
        end
        
    case 'a'
        
        if n == 1 && m == 1
            out.coef = p.coef;
            out.var = p.var;
            out.nvar = p.nvar;
            out.deg = p.deg;
            return
        end
        
        for i = 1:n
            for j = 1:m
                out{i,j}.coef = p(i,j).coef;
                out{i,j}.var = p(i,j).var;
                out{i,j}.nvar = p(i,j).nvar;
                out{i,j}.deg = p(i,j).deg;
            end
        end
end

end


