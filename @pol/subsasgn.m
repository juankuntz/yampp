function out = subsasgn(p,s,q)

% Overloaded subasgn so that it functions as usual, but if we are assigning
% to a matrix, that is using '()', the function makes sure that all of the
% entries of the assigned matrix are written in all relevant independent
% variables. Only matrices supported, no higher dimensional arrays. Also,
% it updates the degrees.

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
        
        if numel(s.subs)>2
            error('Only matrix indexing permitted for polynomials, no higher dimensional arrays allowed. For example, X(3,4) = x(1), is fine, but X(3,4,1) = x(1) is not.');
        end
        
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
            
            % Update the degree.
            
            DEG = max(p(1,1).deg,q(1,1).deg);

            for i = 1:numel(p(:,1))
                for j = 1:numel(p(1,:))
                    p(i,j).deg = DEG;
                    q(i,j).deg = DEG;
                end
            end


        end

        % If p is empty, we need to initialise it as a polynomial and make
        % sure all its entries are written in the same variables (those of
        % q). Also load the degree of q into that of p.
        
        if isempty(p)
            
            clear p
            
            n = s.subs{1};
            m = 1;
            if numel(s.subs) > 1
                m = s.subs{2};
            end
            
            p = pol(zeros([n,m]));
            if ~isempty(q.var.symb)     % If p has variables, load them into q.
                for i = 1:n
                    for j = 1:m
                        p(i,j).var = q.var;
                        p(i,j).deg = q.deg; % Update degree.
                    end
                end
            end
        end
        
        % If q is a double convert it into a polynomial in the variables of
        % p.
        
        if isa(q,'double') && ~isempty(q)
            temp = pol(q); 
            clear q;
            q = temp;
            if ~isempty(p(1,1).var.symb)
                for i = 1:numel(q(:,1))
                    for j = 1:numel(q(1,:))
                        q(i,j).var = p.var; 
                        q(i,j).deg = p.deg; % Update degree.
                    end
                end
            end
        end
        
        % Now we can call the builtin subsasgn.
        
        out = builtin('subsasgn',p,s,q);    
        
        % If q is empty, we may have to update the list of variables, thus
        % run cleanpol.
        
        if isempty(q)
            clear temp;
            temp = out; 
            clear out;
            out = cleanpol(temp);
        end
    case '{}'
        out = builtin('subsasgn',p,s,q);
    case '.'
        out = builtin('subsasgn',p,s,q);
end
end
