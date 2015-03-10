function vars(s,n) 

% Creates in the workspace of the caller of vars the vector independent
% variables [s(1),s(2),...,s(n)]^T. The function will return an error if
% variables s are already declared. The symbol of the variable must be 
% a single character.

% Alternatively, one can call the function by typing
% >> vars x n
% instead of vars('x',n);

% Juan Kuntz, 10/03/2015

temp = pol; % Declare zero polynomial.

if ~isa(s,'char') || ~isscalar(s)
    error('Error: The symbol s of the independent variable must be a single character.')
end 

% Check whether variable s already exists in the workspace of the function 
% in which this new variables are being declared, if this is the case, 
% abort and return an error.

if evalin('caller',['exist(',sprintf( '\''' ),s,sprintf( '\''' ),')'])
    error(['Error: There already is some variable ',s,' defined in the workspace of the function that is trying to define these new variables.']);
end

% If the number of variables is specified by a string convert it into a
% double.

if isa(n,'char')
    temp = n; clear n; n = str2double(temp); clear temp
end

% Check that n is an integer.

if n-floor(n) ~= 0
    error('Error: The number n of independent variables must be an integer.')
end

obj(n,1) = pol; % Declare vector of n zero polynomials.

% Turn obj into vector of variables.

tempchoose = ncktab(n+1);
for i = 1:n
    temp.symb = s; temp.ncomp = n;
    obj(i,1).var = temp;
    obj(i,1).coef = [1;n+2-i]; % Declare coefficients in grlex order. However, this must be done after the variables are declared (otherwise the set.var function jumbles the coefficients).
    obj(i,1).deg = 1;
    obj(i,1).nvar = n;
    obj(i,1).choose = tempchoose;
    clear temp
end
clear tempchoose;

% Save obj as variable s in workspace of function declaring these
% variables.

assignin('caller',s,obj);

end