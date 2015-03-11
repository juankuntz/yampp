function out = subsasgn(p,s,q)

% Overloaded subasgn so that it functions as usual, but if we are assigning
% to a matrix, that is using '()', the function makes sure that all of the
% entries of the assigned matrix are written in all relevant independent
% variables.

% For example,
%
% >> vars x 3
% >> u = x*x';
% >> vars y 1
% >> u(1,1) = y;
% >> u(1,1).var
%
% returns
% 
% u(1,1).var.symb = 'xy'
% u(1,1).var.ncomp = [3,1]
%
% instead of 
%
% u(1,1).var.symb = 'x'
% u(1,1).var.ncomp = 3
%
% Juan Kuntz, 11/03/2015.

switch s.type
    
    case '()'
        
        % If q and p are both polynomials, then make sure they have the
        % same variables.
            
        if isa(q,'pol') && isa(p,'pol') && ~strcmp(p(1,1).var.symb,q(1,1).var.symb)
            [pnotq,qnotp] = varcomp(p(1,1),q(1,1));  

            if ~isempty(pnotq.symb)
                for i = 1:numel(q(:,1))
                    for j = 1:numel(q(1,:))
                        q(i,j).var = pnotq; 
                    end
                end
            end
            
            if ~isempty(qnotp.symb)
                for i = 1:numel(p(:,1))
                    for j = 1:numel(p(1,:))
                        p(i,j).var = qnotp;
                    end
                end
            end

        end

        % If p is empty, we need to initialise it as a polynomial in order
        % to use the buildin subsasgn below.
        
        if isempty(p)
           p = pol; 
           p.var = q(1,1).var;
        end
        
        % If q is a double convert it into a polynomial in the variables of
        % p.
        
        if isa(q,'double')
            temp = pol(q); 
            clear q;
            q = temp;
            for i = 1:numel(q(:,1))
                for j = 1:numel(q(1,:))
                    q(i,j).var = p.var; 
                end
            end
        end
        
        % Now we can call the builtin subsasgn.
        
        out = builtin('subsasgn',p,s,q);    
    case '{}'
        out = builtin('subsasgn',p,s,q);
    case '.'
        out = builtin('subsasgn',p,s,q);
end
end
