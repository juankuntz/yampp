classdef scp < matlab.mixin.SetGet
    
    % Class definition for scp: meant to model moment problems (scp is 
    % short for "sequence conic program").
    
    % Juan Kuntz, 16/02/2015, last edited 15/05/2015.
    
    properties 

        % Properties describing the moment problem. Remark: we only allow 
        % the user to add extra constraints, not delete or replace 
        % previously declared constraints. Similarly for the polynomial
        % equalities and inequalities describing the support of the
        % measure.
        
        nvar = [];          % Dimension of underlying space.
        obj = {[],[]};      % Contains the description of the objective/s.
        
        supineq = [];       % Polynomial ineqequalities defining the support of the measure.
        supeq = []          % Polynomial equalities defining the support of the measure -- not implemented yet.
        
        eqcon = {[],[]};         % Linear equality constraints on the sequence of moments of the measure.
        ineqcon = {[],[]};       % Linear inequality constraints on the sequence of moments of the measure -- not implemented yet.
        
        % Properties describing the relaxations to be solved.
        
        rtyp = [];      % Types of relaxation {D,DD,SDD,{FKW,k},PSD} to be solved.
        rord = [];      % Orders of relaxation to be solved.
        rlst = [];      % Cell array containing the specification of each 
                        % relaxation to be solved, used in combination with 
                        % mkscp by solvescp to build the relaxations.
        mult = pol(1);       % Multiplier s
        multpow = @(d) 1; % Multiplier power function d |-> multpow(d), so that the multiplier in relaxations of order d is s^multpow(d).
        
        % Solver options
        
        ops = []; % Solver options.
        parallelise = 0; % Parallelise option.
        
        % Property storing data on solved relaxations
        
        sol = []; % A cell array, each entry of which contains all the data on one solved relaxation.

        % This property requires fixing.
        
        dualres = 0;
        
    end
    
    methods
        
        % Class constructor.
        
        function obj = scp
            
            % Set default options to sedumi with default yalmip options.
            
            set(obj,'ops',sdpsettings('solver','sedumi'));

        end
        
        % Property set methods with error messages. Remark: we only allow 
        % the user to add extra constraints, not delete or replace 
        % previously declared constraints. Similarly for the polynomial
        % equalities and inequalities describing the support of the
        % measure.
        
        function set.ops(cp,data)
            try 
                
            if isstruct(data)
                cp.ops = {data,data,data,data,data,data};
            else
                ops = []; label =[];
                for i = 1:numel(data)
                    if ischar(data{i})
                       label{end+1} = data{i};
                    else
                       if ~isempty(ops)
                            error
                        end
                       ops = data{i};
                    end
                end
                
                if isempty(ops) || isempty(label)
                        error
                end
                    
                for i = 1:numel(label)
                        temp = cp.ops;
                    switch label{i}
                        case {'d','D'}
                            cp.ops = {ops,temp{2,end},temp{3,end},temp{4,end},temp{5,end},temp{6,end}};
                        case {'dd','DD'}
                            cp.ops = {temp{1,end},ops,temp{3,end},temp{4,end},temp{5,end},temp{6,end}};
                        case {'sdd','SDD'}
                            cp.ops = {temp{1,end},temp{2,end},ops,temp{4,end},temp{5,end},temp{6,end}};
                        case {'fwk','FWK'}
                            cp.ops = {temp{1,end},temp{2,end},temp{3,end},ops,temp{5,end},temp{6,end}};
                        case {'psd','PSD'}
                            cp.ops = {temp{1,end},temp{2,end},temp{3,end},temp{4,end},ops,temp{6,end}};
                        case {'nn','NN'}
                            cp.ops = {temp{1,end},temp{2,end},temp{3,end},temp{4,end},temp{5,end},ops};
                    end
                end
                    
            end
            
            catch
                error('Invalid format of input specifying the solver options for the various relaxations.');
            end
                
            
        end
        
        function set.rtyp(cp,data)
            try
                types = cp.rtyp;
                for i = 1:numel(data)
                   s = data{i};
                   if ischar(s) && (strcmpi(s,'D') || strcmpi(s,'DD') || strcmpi(s,'SDD') || strcmpi(s,'PSD') || strcmpi(s,'NN'))
                       flg = 1;
                       for j = 1:numel(types)
                           if ischar(types{j}) && istrcmp(s,types{j})
                              flg = 0; 
                           end
                       end
                   elseif ischar(s{1}) && strcmpi(s{1},'FWK') && numel(s) == 2
                       flg = 1;
                       for j = 1:numel(types)
                           type = types{j};
                           if ~ischar(type) && s{2} == type{2};
                              flg = 0; 
                           end
                           clear type
                       end
                   else
                       error;
                   end
                   
                   if flg
                       cp.rtyp{end+1} =  s;
                   end
                   
                   clear s
                end
                
            catch
                error('Invalid format of input specifying the type of relaxations to be sovled.');
            end
        end
        
        function set.obj(cp,data)
            
            try 
                % Declare shorthands.
                
                c = data{1};
                f = data{2};
                
                % Check that data is of the correct format.
                
                flg = 0;
                for i=1:numel(c(:,1))
                    if ~isa(c(i,:),'char') || (strcmp('inf',c(i,:)) && strcmp('sup',c(i,:)))
                        flg = 1;
                    end
                end
                
                % If the objective function is actually specified as a 
                % double, then convert it into a pol (a user may do this if
                % she is optimising the mass of a meaure.
                
                if isa(f,'double')
                    f = pol(f);
                end
                
                if ~isa(f,'pol') || ~isequal(numel(f),numel(c(:,1))) || numel(data) > 2 || flg
                    error;
                end
                
                % If everything is fine, store the constraints.
                
                cp.obj{1} = [cp.obj{1};c];
                cp.obj{2} = [cp.obj{2};f(:)];
                
            catch % If data has been specified poorly, return error message.
                error('The objective must be specified by a cell C containing two objects, a matrix chars c = C{1} and vector, with its length being the number of rows of c, of pols p = C{2}. Each row of c must contain either inf or sup. This adds the objectives c(i,:) p(i) (elementwise) to the program.');
            end   
        end 
        
        
        function set.supeq(cp,p)
            
            if ~isa(p,'pol') || min(size(p)) ~= 1 || numel(data) > 1
                error('The support of the measure must be specified by vectors of polynomials.');
            end
            
            cp.supeq = [cp.supeq,p(:)];
        end
        
        function set.supineq(cp,p)
            
                if ~isa(p,'pol') || min(size(p)) ~= 1 
                    error('The support of the measure must be specified by vectors of polynomials.');
                end
                
                cp.supineq = [cp.supineq,p(:)];
        end 
        
        function set.eqcon(cp,data)
            
            try 
                % Declare shorthands.
                
                if isa(data,'pol')
                    p = data;
                    c = zeros(size(p));
                elseif numel(data) == 1 
                    p = data{1};
                    c = zeros(size(p));
                else
                    p = data{1};
                    c = data{2};
                end
                
                % If the polynomial defining the constraint is actually 
                % specified as a double, then convert it into a pol (a user
                % may do this if she is constraining the mass of a meaure).
                
                if isa(p,'double')
                    p = pol(p);
                end
                
                % Check that data is of the correct format.
                
                if ~isa(p,'pol') || ~isa(c,'double') || ~isequal(size(p),size(c)) || (~isa(data,'pol') && numel(data) > 2)
                    error;
                end
                
                % If everything is fine, store the constraints.
                
                cp.eqcon{1} = [cp.eqcon{1};p(:)];
                cp.eqcon{2} = [cp.eqcon{2};c(:)];
                
            catch % If data has been specified poorly, return error message.
                error('The equality constraints must be specified by a cell C containing two objects, a matrix of polynomials p = C{1} and a matrix, of the same dimensions, of doubles c = C{2}. This adds the constraint integral(p) == c (elementwise) to the program. If c = 0, c may be ommitted from C.');
            end   
        end 
            
        function set.ineqcon(cp,data)
            
            % Remark: we only allow the user to add extra constraints, not
            % delete or replace previously declared constraints.
            
            try 
                % Declare shorthands.
                
                if isa(data,'pol')
                    p = data;
                    c = zeros(size(p));
                elseif numel(data) == 1
                    p = data{1};
                    c = 0;
                else
                    p = data{1};
                    c = data{2};
                end
                
                % If the polynomial defining the constraint is actually 
                % specified as a double, then convert it into a pol (a user
                % may do this if she is constraining the mass of a meaure).
                
                if isa(p,'double')
                    p = pol(p);
                end
                
                % Check that data is of the correct format.
                
                if ~isa(p,'pol') || ~isa(c,'double') || ~isequal(size(p),size(c)) || (~isa(data,'pol') && numel(data) > 2)
                    error;
                end
                
                % If everything is fine, store the constraints.
                
                cp.ineqcon{1} = [cp.ineqcon{1};p(:)];
                cp.ineqcon{2} = [cp.ineqcon{2};c(:)];
                
            catch % If data has been specified poorly, return error message.
                error('The inequality constraints must be specified by a cell C containing two objects, a matrix of polynomials p = C{1} and a matrix, of the same dimensions, of doubles c = C{2}. This adds the constraint integral(p) >= c (elementwise) to the program. If c = 0, c may be ommitted from C.');
            end   
        end 
        
        % Property get methods.
        
        function out = get.ops(cp)

            out = {'D',cp.ops{1}.solver,cp.ops{1};
                   'DD',cp.ops{2}.solver,cp.ops{2};
                   'SDD',cp.ops{3}.solver,cp.ops{3};
                   'FWK',cp.ops{4}.solver,cp.ops{4};
                   'PSD',cp.ops{5}.solver,cp.ops{5};
                   'NN',cp.ops{6}.solver,cp.ops{6}};
        end
    end
end
