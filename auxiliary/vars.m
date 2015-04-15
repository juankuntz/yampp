function vars(s,n) 

% Creates in the workspace of the caller of vars the vector independent
% variables [s(1),s(2),...,s(n)]^T. The function will return an error if
% variables s are already declared. The symbol of the variable must be 
% a single character.

% Alternatively, one can call vars by typing
%
% >> vars x n
%
% instead of 
%
% >> vars('x',n);

% Juan Kuntz, 10/03/2015

temp = pol; % Declare zero polynomial.

if ~isa(s,'char') || ~isscalar(s)
    error('The symbol s of the independent variable must be a single character.')
end 

% Check whether variable s already exists in the workspace of the function 
% in which this new variables are being declared, if this is the case, 
% abort and return an error.

if evalin('caller',['exist(',sprintf( '\''' ),s,sprintf( '\''' ),')'])
    error(['There already is some variable ',s,' defined in the workspace of the function that is trying to define these new variables.']);
end

% If the number of variables is specified by a string convert it into a
% double.

if isa(n,'char')
    temp = n; clear n; n = str2double(temp); clear temp
end

% Check that n is an integer.

if n-floor(n) ~= 0
    error('The number n of independent variables must be an integer.')
end

obj(n,1) = pol; % Declare vector of n zero polynomials.

% Turn obj into vector of variables, since the properties of pols are 
% private, the code to do this needs to be stored in a method of pol; 
% devar.m.

out = decvar(obj,s,n);

% Save obj as variable s in workspace of function declaring these
% variables.

assignin('caller',s,out);

end