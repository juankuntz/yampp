function solvescp(cp)

%
% Juan Kuntz, 16/02/2015, last edited 14/03/2015.

% Extract dual point and save in sol.
% 
% Compute, by hand (as in http://users.isy.liu.se/en/rt/johanl/2009_OMS_DUALIZE.pdf), the primal and dual residues, use dual() to extract dual variables.
% 
% Compare these with the residues returned by check().

% Missing error/warning messages.


% Check that the degree of no polynomial in the problem's description
% exceeds (2-times) the relaxation order. 

if (~isempty(cp.supineq) && 2*min(cp.rord) < prop(cp.supineq(1),'d')) || (~isempty(cp.eqcon{1}) && 2*min(cp.rord) < prop(cp.eqcon{1}(1),'d')) || (~isempty(cp.ineqcon{1}) && 2*min(cp.rord) < prop(cp.ineqcon{1}(1),'d')) || (~isempty(cp.obj{2}) && 2*min(cp.rord) < prop(cp.obj{2}(1),'d')) 
    error('The degree of one of the polynomial defining the support is bigger than the order of the moments included in the moment problem, increase the relaxation order of the problem.'); 
end

setvars(cp);    % Decide the dimension of the underlying space.
mklst(cp);      % Make list containing the specifications of each relaxations to be solved.

% Build scp.

sol{numel(cp.rlst)} = [];

if cp.parallelise == 1
    parfor i = 1:numel(cp.rlst)
        sol{i} = solve(cp,mkrel(cp,cp.rlst{i}));
    end
else
    for i = 1:numel(cp.rlst)
        sol{i} = solve(cp,mkrel(cp,cp.rlst{i}));
    end
end

for i = 1:numel(cp.rlst)
   cp.sol{i} = sol{i}; 
end

end


function sol = solve(cp,rel)

sol = rel; 
clear rel;

% Shorthands

n = cp.nvar; d = sol.relorder;

% Fish out solver options.

switch sol.reltype
    case {'d','D'}
        ops = cp.ops{1,end};
    case {'dd','DD'}
        ops = cp.ops{2,end};
    case {'sdd','SDD'}
        ops = cp.ops{3,end};
    case {'fwk','FWK'}
        ops = cp.ops{4,end};
    case {'psd','PSD'}
        ops = cp.ops{5,end};
end

if strcmp(sol.minmax,'inf')
    temp = optimize(sol.ycons,sol.yobj,ops); % Compute solution.

    sol.dval = -value(dual(sol.ycons(1))); % Not sure why we need the minus sign here...
elseif strcmp(sol.minmax,'sup')
    temp = optimize(sol.ycons,-sol.yobj,ops);

    sol.dval = value(dual(sol.ycons(1))); 
end

    % Store solution.

    sol.pval = value(sol.yobj);
    sol.ppoint = seq(n,2*d,value(sol.yvar));

    % Store primal residues in a formated table.

    [temp2,~] = check(sol.ycons);

    sol.pres{1,2} = [];

    l = 0;
    if ~isempty(cp.mass)
        sol.pres{1,1} = 'Mass constraint.';
        sol.pres{1,2} = temp2(1);
        l = l + 1;
    end

    if ~isempty(cp.eqcon{1})
        if l == 0
            sol.pres{1,1} = 'Equality constraints.';
            sol.pres{1,2} = temp2(1);
        else
            sol.pres{2,1} = 'Equality constraints.';
            sol.pres{2,2} = temp2(2);
        end
        l = l + 1;
    end

    if l == 0
        sol.pres{1,1} = 'Cone constraints.';
        sol.pres{1,2} = temp2(1);
    else
        sol.pres{end+1,1} = 'Cone constraints.';
        sol.pres{end,2} = temp2(l+1);
    end
    l = l + 1;

   for j = 1:numel(cp.supineq)
        sol.pres{end,2} = [sol.pres{end,2};temp2(l+j)];
   end

    % Store yalmip output info.

    sol.info = temp; clear temp temp2;

    % If problem is unbounded or infeasible store the appropiate +-inf
    % into the primal value.

    if sol.info.problem == 1 
        if strcmp(sol.minmax,'inf')
            sol.pval = inf;
        else
            sol.pval = -inf;
        end
    elseif sol.info.problem == 2
        if strcmp(sol.minmax,'inf')
            sol.pval = -inf;
        else
            sol.pval = inf;
        end
    end

    % If dual residues not requested, exit (saves considerable time).

    if cp.dualres == 0
        return
    end

    % Compute and store the dual residues in a formated table.      

    if ~isempty(cp.mass) && ~isempty(cp.eqcon{1})
        t = dual(sol.ycons(1)); % Dual of the mass constraint.
        t(2:numel(cp.eqcon{1})+1,1) = dual(sol.ycons(2)); % Dual of the equality constraints
    elseif ~isempty(cp.mass)
        t = dual(sol.ycons(1)); % Dual of the mass constraint.
    elseif ~isempty(cp.eqcon{1})
        t = dual(sol.ycons(1)); % Dual of the equality constraints
    end

    temp = sol.F'*t; 
    clear t;


    if strcmp(sol.minmax,'inf')
        temp = temp + sol.b;
    else
        temp = temp - sol.b;
    end

    % Onto the support constraints.

    sol.dres{1,1} = 'Equality constraints';

    switch sol.reltype
        case {'D','DD'}

            for j = 1:numel(cp.supcon)+1
                if isempty(cp.mass) && isempty(cp.eqcon{1})
                    X{j} = dual(sol.ycons(j));
                elseif isempty(cp.mass) || isempty(cp.eqcon{1})
                    X{j} = dual(sol.ycons(1+j));
                else
                    X{j} = dual(sol.ycons(2+j));
                end
            end

            for j = 1:numel(cp.supcon)+1
                temp = temp + sol.A{j}*X{j};
            end

            sol.dres{1,2} = -max(abs(temp));   % Return Linfty norm of equality constraing violations
            clear temp

            sol.dres{2,1} = 'Cone constraints';
            for j = 1:numel(cp.supcon)+1
                sol.dres{2,2} = [sol.dres{2,2};min(X{j})];
            end

        case 'SDD'

            if isempty(cp.mass) && isempty(cp.eqcon{1})
                L = 0;
            elseif isempty(cp.mass) || isempty(cp.eqcon{1})
                L = 1;
            else
                L = 2;
            end

            % Initialise SOC dual variables

            X{1} = zeros(3*nchoosek(n+d,d),1);
            
            for j = 1:numel(cp.supcon)
                dc = deg(cp.supcon(j));
                X{j+1} = zeros(3*nchoosek(n+floor(d-dc/2),floor(d-dc/2),1));
            end

            % Populate the SOC duals

            J = 1;
            for j = 1:numel(cp.supcon)+1
                temp2 = sol.A{j};
                for k = 1:numel(temp2)
                    X{j}(1+3*(k-1):3*k,:) = dual(sol.ycons(L+J));
                    temp = temp + temp2{k}*X{j}(1+3*(k-1):3*k,:);
                    J = J + 1;
                end
                clear temp2
            end

            sol.dres{1,2} = -max(abs(temp));    % Return Linfty norm of equality constraing violations
            clear temp

            sol.dres{2,1} = 'Cone constraints';
            sol.dres{2,2} = [];
            for j = 1:numel(cp.supcon)+1
                temp = X{j};
                sol.dres{2,2} = [sol.dres{2,2}; X{j}(1)-sqrt(X{j}(2)^2+X{j}(3)^2)];
                for k = 2:numel(X{j})/3
                    sol.dres{2,2} = min(sol.dres{2,2},X{j}((k-1)*3+1)-sqrt(X{j}((k-1)*3+2)^2+X{j}((k-1)*3+3)^2));
                end
            end

        case 'FKW'

        case 'PSD'

            for j = 1:numel(cp.supcon)+1
                if isempty(cp.mass) && isempty(cp.eqcon)
                    X{j} = dual(cp.ycons(j));
                elseif isempty(cp.mass) || isempty(cp.eqcon)
                    X{j} = dual(sol.ycons(1+j));
                else
                    X{j} = dual(sol.ycons(2+j));
                end
            end

            for j = 1:numel(cp.supcon)+1
                temp2 = [];
                A = sol.A{j};
                for k = 1:numel(sol.A{j})
                    temp2 = [temp2;sum(sum(X{j}.*A{k}))];
                end
                clear A;

                temp = temp + temp2; clear temp2
            end

            sol.dres{1,2} = -max(abs(temp)); % Return Linfty norm of equality constraing violations
            clear temp

            sol.dres{2,2} = [];

            for j = 1:numel(cp.supcon)+1
                sol.dres{2,1} = 'Cone constraints';
                sol.dres{2,2} = [sol.dres{2,2};min(eig(X{j}))];
            end
    end
end
