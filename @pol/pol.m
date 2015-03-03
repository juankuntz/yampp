classdef pol
    % WE SHOULD ADD SOMETHING THAT CHECKS THAT isymbol DOES NOT HAVE
    % REPEATED SYMBOLS. ALSO THAT THE BASIS MATRIX ONLY CONTAINS INTEGERS,
    % ETC. ALSO NO REPEATED MONOMIALS IN THE BASIS.
    
    % ALSO SAME NUMBER OF COEFFICIENTS AS OF BASIS.
    
    % MAYBE ALSO ADD SOMETHING THAT AUTOMATICALLY GETS RID OF ZERO
    % COEFFICIENTS AND THE CORRESPONDING ELEMENTS OF THE BASIS (THAT WAY WE
    % CAN REDUCE THE MEMORY EACH OF THESE THINGS TAKE UP). WATCH OUT THIS
    % WILL SCREW WITH combinepols.m.
    
    % ALSO EACH ELEMENT IN AN ARRAY MUST HAVE SAME BASIS, SYMBOL ETC?

    properties %(SetAccess = private)
        coef = []; % 2 x length(coef) array. The first row contains a nonzero coefficient of the polynomial and the corresponding entry of the second row contains the grlex (in sum(var.ncomp) variables) rank.
        var = [];  % Structure that contains two arrays, symb and ncomp; symb contains the symbols of the independent variables (single characters).
        nvar = []; % Total number of variables; so we don't need to compute it repeatedly.
        deg = []; % Total degree; so we don't need to compute it repeatedly.
        choose = []; % (i,j) contains i choose j; so we do not have to call nchoosek repeatedly (it's expensive).
    end 
    
    methods % Add method, find basis element -- returns where in the basis said element is.
            % Add method var2loc    % Given isymbol and number of the symbol, returns what row in the basis matrix corresponds to said variable.
            
            function obj = pol(varargin) % Constructor method
                
             % If only one argument is specified, construct constant
             % polynomial. If no argument is specified, construct 0
             % polynomial. This is used internally and to convert doubles
             % to pols.
             
            if nargin <= 1;
                obj.deg = 0;
                obj.nvar = 0;
                temp.symb = [];
                temp.ncomp = [];
                obj.var = temp;
                clear temp
                if nargin == 1  % Zero polynomial is identified with the one that has empty coef array.
                    obj.coef = [varargin{1};1];
                end
                return
            end
            
            % Otherwise we are declaring variables.
            
            s = varargin{1}; n = varargin{2};
            
            if ~isa(s,'char') || ~isscalar(s)
                disp('Error: The symbol s of the independent variable must be a single character.')
                return
            end 
            if isa(n,'char')
                temp = n; clear n; n = str2double(temp); clear temp
            end
            if n-floor(n) ~= 0
                disp('Error: The number c of independent variables must be an integer.')
                return
            end
            obj(n,1) = pol;
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
            
            % MAKE SURE VARIABLE S DOES NOT EXIST ALREADY, USE EVALIN
            % FUNCTION FOR THIS.
            
            if evalin('caller',['exist(',sprintf( '\''' ),s,sprintf( '\''' ),')'])
                disp(['Error: There already is some variable ',s,' defined in the workspace of the function that is trying to define these new variables.']);
                return
            end
            assignin('caller',s,obj);
            obj = pol;
            end
        
            function obj = set.var(obj,new)
               
                if isempty(obj.var) % If there are no variables already, then the variables are just the new ones.
                    [obj.var.symb,I] = sort(new.symb); % Order the variables alphabetically.
                    obj.var.ncomp = new.ncomp(I);
                    return
                end
                
                % Otherwise mix old and new variables.
                oldsymb = obj.var.symb;
                tempsymb = [oldsymb,new.symb]; 
                nsymbnew = numel(tempsymb); % Need this later
                obj.var.symb = []; 
                
                tempncomp = [obj.var.ncomp,new.ncomp]; 
                nvarsold = sum(obj.var.ncomp); % Need this later
                obj.var.ncomp = [];


                [obj.var.symb,I] = sort(tempsymb); % Order the variables alphabetically.
                obj.var.ncomp = tempncomp(I);
                clear tempncomp tempsymb
                
                [obj,oldchoose] = updatenvar(obj,sum(obj.var.ncomp)); % Update nchoose and nvar properties, also store old choose table for next step.
                obj = updatecoef(obj,nvarsold,nsymbnew,oldsymb,oldchoose); % Rewrite coefficients in new ranking.
            end
    end
end
