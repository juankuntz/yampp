function out = subsref(p,s)

% Overloaded so that with bracket index cleanpol is called in case we've
% extracted polynomials in less variables or of lower power. However, we
% can no longer access the properties of a polynomial directly, we will
% need to implement other methods to do this. This is because MATLAB's
% special way of deciding how many outputs this function should have which
% eludes me (and it appears the MATLAB documentation too). Some relevant
% info See http://www.mathworks.com/matlabcentral/answers/101955-why-do-i-receive-errors-when-overloading-subsref-for-types-and-for-matlab-classes

% Juan Kuntz, 12/03/2015.

% Remeber: Inside the methods of a class, the short hands for subsref and
% subasgn always call the build in subsref and subsasgn, see documentation.
% To call this one we need to use the explicit call subsref(A,s).

switch s.type
    
    case '()'
        if numel(s.subs)>2
            error('Only matrix indexing permitted for polynomials, no higher dimensional arrays allowed. For example, X(3,4) = x(1), is fine, but X(3,4,1) = x(1) is not.');
        end
        
        out = builtin('subsref',p,s);  
        out = cleanpol(out);

    case '{}'
        % Not implemented.
    case '.'
        % Not implemented.
end

end