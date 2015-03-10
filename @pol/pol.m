classdef pol
    
    % The class definition of object pol used to represent polynomials.
    % Essentially it represents it a polynomial by storing its nonzero
    % coefficients and their rank in the graded lexigraphic ordering (grlex)
    % into the property coef. It also saves the variables of the
    % polynomial.
    
    % Juan Kuntz, 05/02/2015, last edited 10/03/2015

    properties (SetAccess = private)
        coef = []; % 2 x length(coef) array. The first row contains a nonzero coefficient of the polynomial and the corresponding entry of the second row contains the grlex (in sum(var.ncomp) variables) rank.
        var = [];  % Structure that contains two arrays, symb and ncomp; symb contains the symbols of the independent variables (single characters).
        nvar = []; % Total number of variables; so we don't need to compute it repeatedly.
        deg = []; % Total degree; so we don't need to compute it repeatedly.
    end 
    
    methods             
            function obj = pol(varargin) 
             
            % Constructor method: can be used to declare either zero
            % polynomials, if no argument is specified, or constant
            % polynomials if an argument is specified

            if nargin > 1 || (nargin == 1 && ~isdouble(varargin{1}))
                error('The constructor method of pol objects either takes in no arguments, in which case it returns the zero polynomial, or it takes a single double that is used to declare the constant polynomial with that double as its coefficient');
            end
            
            obj.deg = 0;
            obj.nvar = 0;
            temp.symb = [];
            temp.ncomp = [];
            obj.var = temp;
            clear temp
            if nargin == 1  % Zero polynomial is identified with the one that has empty coef array.
                obj.coef = [varargin{1};1];
            end

            end
        
            function obj = set.var(obj,new)
                
                % The set method for var is used to add extra variables. It
                % is not possible to remove or replace variables.
               
                if isempty(obj.var) % If there are no variables already, then the variables are just the new ones.
                    [obj.var.symb,I] = sort(new.symb); % Order the variables alphabetically.
                    obj.var.ncomp = new.ncomp(I);
                    return
                end
                
                % Otherwise mix old and new variables.
                
                oldsymb = obj.var.symb;
                tempsymb = [oldsymb,new.symb]; 
                nsymbnew = numel(tempsymb); % Need this later when updating the coefficients.
                obj.var.symb = []; 
                
                tempncomp = [obj.var.ncomp,new.ncomp]; 
                nvarsold = sum(obj.var.ncomp); % Need this later when updating the coefficients.
                obj.var.ncomp = [];

                [obj.var.symb,I] = sort(tempsymb); % Order the variables alphabetically.
                obj.var.ncomp = tempncomp(I);
                clear tempncomp tempsymb
                
                obj = updatenvar(obj,sum(obj.var.ncomp)); % Update nvar property, also store old choose table for next step.
                
                obj = updatecoef(obj,nvarsold,nsymbnew,oldsymb); % Rewrite coefficients in new ranking.
            end
    end
end
