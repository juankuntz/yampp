function solvescp(cp)

% Comments below are bogus, only single argument cp that has already been
% set.
%
% Function called by user to solve the moment problem. Three arguments must
% be specified:
%
% 1 - The moment problem (short for conic program).
% 2 - The order of the relaxation (overload this to allow for several orders?).
% 3 - A string specifying the type of relaxation to be solved. This can be
% one of the following:
%
% -- D (not implemented yet).
% -- DD (not implemented yet).
% -- SDD (not implemented yet).
% -- FWK (not implemented yet). If chosen must specify factor width in
% fourth argument.
% -- PSD: The traditional moment conditions (i.e. feasible points are
%         contained in the dual of the SOS cone.
% -- ALL (not implemented yet).
%
% Juan Kuntz, 16/02/2015, last edited 16/02/2015.

% Extract dual point and save in sol.
% 
% Compute, by hand (as in http://users.isy.liu.se/en/rt/johanl/2009_OMS_DUALIZE.pdf), the primal and dual residues, use dual() to extract dual variables.
% 
% Compare these with the residues returned by check().

% Missing error/warning messages.

% Build scp.

mkscp(cp);

% Shorthands

n = cp.nvar; d = cp.relorder;

for i = 1:numel(cp.obj)
    if strcmp(cp.minmax,'min')
        temp = optimize(cp.ycons,cp.yobj(i),cp.ops); % Compute solution.
        cp.sol{end+1}.minmax = 'min';
        
        cp.sol{end}.dval = -value(dual(cp.ycons(1))); % Not sure why we need the minus sign here...
    elseif strcmp(cp.minmax,'max')
        temp = optimize(cp.ycons,-cp.yobj(i),cp.ops);
        cp.sol{end+1}.minmax = 'max';
        
        cp.sol{end}.dval = value(dual(cp.ycons(1))); 
    end

        % Store solution.
        
        cp.sol{end}.reltype = cp.reltype;
        cp.sol{end}.obj = cp.obj(i);
        
        cp.sol{end}.ycons = cp.ycons;
        cp.sol{end}.relorder = d;
        
        cp.sol{end}.pval = value(cp.yobj(i));
        cp.sol{end}.ppoint = seq(n,2*d,value(cp.yvar));
        
        
        % Store primal residues in a formated table.
        
        [temp2,~] = check(cp.ycons);
        
        cp.sol{end}.pres{1,2} = [];
        
        l = 0;
        if ~isempty(cp.mass)
            cp.sol{end}.pres{1,1} = 'Mass constraint.';
            cp.sol{end}.pres{1,2} = temp2(1);
            l = l + 1;
        end
        
        if ~isempty(cp.seqeqcon)
            if l == 0
                cp.sol{end}.pres{1,1} = 'Equality constraints.';
                cp.sol{end}.pres{1,2} = temp2(1);
            else
                cp.sol{end}.pres{2,1} = 'Equality constraints.';
                cp.sol{end}.pres{2,2} = temp2(2);
            end
            l = l + 1;
        end
        
        if l == 0
            cp.sol{end}.pres{1,1} = 'Cone constraints.';
            cp.sol{end}.pres{1,2} = temp2(1);
        else
            cp.sol{end}.pres{end+1,1} = 'Cone constraints.';
            cp.sol{end}.pres{end,2} = temp2(l+1);
        end
        l = l + 1;
        
       for j = 1:numel(cp.supcon)
            cp.sol{end}.pres{end,2} = [cp.sol{end}.pres{end,2};temp2(l+j)];
       end
        
        % Store yalmip output info.
        
        cp.sol{end}.info = temp; clear temp;
        
        % If problem is unbounded or infeasible store the appropiate +-inf
        % into the primal value.
        
        if cp.sol{end}.info.problem == 1 
            if strcmp(cp.minmax,'min')
                cp.sol{end}.pval = inf;
            else
                cp.sol{end}.pval = -inf;
            end
        elseif cp.sol{end}.info.problem == 2
            if strcmp(cp.minmax,'min')
                cp.sol{end}.pval = -inf;
            else
                cp.sol{end}.pval = inf;
            end
        end

        % Compute and store the dual residues in a formated table.
        
        if ~isempty(cp.mass) && ~isempty(cp.seqeqcon)
            t = dual(cp.ycons(1)); % Dual of the mass constraint.
            t(2:numel(cp.seqeqcon)+1,1) = dual(cp.ycons(2)); % Dual of the equality constraints
        elseif ~isempty(cp.mass)
            t = dual(cp.ycons(1)); % Dual of the mass constraint.
        elseif ~isempty(cp.seqeqcon)
            t = dual(cp.ycons(1)); % Dual of the equality constraints
        end
        
        temp = cp.F'*t; 
        clear t;

            
        if strcmp(cp.minmax,'min')
            temp = temp + cp.b{i};
        else
            temp = temp - cp.b{i};
        end

        % Onto the support constraints.

        cp.sol{end}.dres{1,1} = 'Equality constraints';
            
    switch cp.reltype
        case {'D','DD'}
            
            for j = 1:numel(cp.supcon)+1
                if isempty(cp.mass) && isempty(cp.seqeqcon)
                    X{j} = dual(cp.ycons(j));
                elseif isempty(cp.mass) || isempty(cp.seqeqcon)
                    X{j} = dual(cp.ycons(1+j));
                else
                    X{j} = dual(cp.ycons(2+j));
                end
            end
        
            for j = 1:numel(cp.supcon)+1
                temp = temp + cp.A{j}*X{j};
            end
            
            cp.sol{end}.dres{1,2} = temp;   clear temp
            
%             cp.sol{end}.dres{2,1} = 'Cone constraints';
%             for j = 1:numel(cp.supcon)+1
%                 cp.sol{end}.dres{2,2} = [cp.sol{end}.dres{2,2};min(X{j})];
%             end

        case 'SDD'
            
            if isempty(cp.mass) && isempty(cp.seqeqcon)
                L = 0;
            elseif isempty(cp.mass) || isempty(cp.seqeqcon)
                L = 1;
            else
                L = 2;
            end
            
            J = 1;
            for j = 1:numel(cp.supcon)+1
                X{j} = [];
                A = cp.A{j};
                for k = 1:numel(A(1,:))
                    X{j} = [X{j};dual(cp.ycons(L+J))];
                    J = J + 1;
                end
            end
            
        case 'FKW'

        case 'PSD'
            
            for j = 1:numel(cp.supcon)+1
                if isempty(cp.mass) && isempty(cp.seqeqcon)
                    X{j} = dual(cp.ycons(j));
                elseif isempty(cp.mass) || isempty(cp.seqeqcon)
                    X{j} = dual(cp.ycons(1+j));
                else
                    X{j} = dual(cp.ycons(2+j));
                end
            end
        
            for j = 1:numel(cp.supcon)+1
                temp2 = [];
                A = cp.A{j};
                for k = 1:numel(cp.A{j})
                    temp2 = [temp2;sum(sum(X{j}.*A{k}))];
                end
                clear A;

                temp = temp + temp2; clear temp2
            end

            cp.sol{end}.dres{1,2} = temp; clear temp
            cp.sol{end}.dres{2,2} = [];
            
            for j = 1:numel(cp.supcon)+1
                cp.sol{end}.dres{2,1} = 'Cone constraints';
                cp.sol{end}.dres{2,2} = [cp.sol{end}.dres{2,2};min(eig(X{j}))];
            end
    end
        
end

end
