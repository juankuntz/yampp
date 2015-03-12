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
            % matrices of constant polynomials if a single matrix of
            % doubles is specified. Eg, pol([1,2]) returns a 1 x 2 vector
            % containing the polynomial 1 in its first entry, and 2 in its
            % second entry. 

            if nargin > 1 || (nargin == 1 && ~isdouble(varargin{1}))
                error('The constructor method of pol objects either takes in no arguments, in which case it returns the zero polynomial, or it takes a matrix of doubles that is used to declare the corresponding matrix of constant polynomials');
            end
            
            if nargin ~= 0
                
            [n,m] = size(varargin{1});
            
                if max([n,m])>1
                    for i = 1:n
                        for j = 1:m
                            obj(i,j) = pol(varargin{1}(i,j));
                        end
                    end
                    return
                end
            end
            
            % If we've made it to here, we're now either in the 0
            % polynomial or the scalar constant polynomial case.
            
            obj.deg = 0;
            obj.nvar = 0;
            temp.symb = [];
            temp.ncomp = [];
            obj.var = temp;
            clear temp
                
            % Zero polynomial is identified with the one that has empty
            % coef array, in which case we are done.
            
            if nargin == 0 
                return                
            end

            obj.coef = [varargin{1};1];
            
            end
        
            function obj = set.var(obj,new)
                
                % The set method for var is used to add extra variables. It
                % is not possible to remove or replace variables. The
                % variables in new must be already ordered alphabetically.
               
                if isempty(obj.var) || isempty(obj.var.symb) % If there are no variables already, then the variables obj are just the new ones.
                    obj.var = new;
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
