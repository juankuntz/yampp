function solvescp(cp)

% Juan Kuntz, 16/02/2015, last edited 16/04/2015.

% Extract dual point and save in sol.
% 
% Compute, by hand (as in http://users.isy.liu.se/en/rt/johanl/2009_OMS_DUALIZE.pdf), the primal and dual residues, use dual() to extract dual variables.
% 
% Compare these with the residues returned by check().

% Missing error/warning messages.


% Check that the degree of no polynomial in the problem's description
% exceeds (2-times) the relaxation order. 

if (~isempty(cp.supineq) && 2*min(cp.rord) < prop(cp.supineq(1),'d')) || (~isempty(cp.obj{2}) && 2*min(cp.rord) < prop(cp.obj{2}(1),'d')) 
    error('The degree of one of the polynomials defining the measure`s support or defining the objective is bigger than the order of the moments included in the moment problem, increase the relaxation order.'); 
end

% Check whether any constraints are ignored.

if (~isempty(cp.eqcon{1}) && prop(cp.eqcon{1}(1),'d') > 2*min(cp.rord)) || (~isempty(cp.ineqcon{1}) && prop(cp.ineqcon{1}(1),'d') > 2*min(cp.rord))
    warning('The polynomial specifying at least one equality or inequality constraint is of higher degree than (two times) the order of at least one relaxation. When solving a relaxation of order d, any equality or inequality constraint of degree greater than 2d is ignored.');
end
clear temp1 temp2
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

cp.sol = sol; 
end

function sol = solve(cp,rel)

sol = rel; 
clear rel;

% Shorthands.

n = cp.nvar; d = sol.rord;

% Solve relaxation.

if strcmp(sol.objs,'inf')
    temp = optimize(sol.ycons,sol.yobj,sol.ops); % Compute solution.
elseif strcmp(sol.objs,'sup')
    temp = optimize(sol.ycons,-sol.yobj,sol.ops);
end

% Store primal and dual solutions

sol.pval = value(sol.yobj);
sol.ppoint = seq(n,2*d,value(sol.yvar));
sol.lambda = value(dual(sol.ycons(1)));

if ~isempty(sol.F)
    sol.lambda = value(dual(sol.ycons(1)));
    if strcmp(sol.objs,'inf')
        sol.dval = -sol.f'*value(dual(sol.ycons(1))); % Minus sign here because yalmip fits this to a maximisation problem, see Lofberg's paper on Dualize.
    elseif strcmp(sol.objs,'sup')
        sol.dval = sol.f'*sol.lambda; 
    end
    for i = 1:numel(sol.A)
        sol.X{i} = value(dual(sol.ycons(i+1)));
    end
else 
    sol.dval = 0;
    for i = 1:numel(sol.A)
        sol.X{i} = value(dual(sol.ycons(i)));
    end
end

% Store primal residues in a formated table.

[temp2,~] = check(sol.ycons);

sol.pres{1,2} = [];
l = 0;
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

% Compute and solve dual residues.

if strcmpi(sol.rtyp,'PSD')
    computedualsres(sol);
end

% Store yalmip output info.

sol.info = temp; clear temp temp2;

% If problem is unbounded or infeasible store the appropiate +-inf
% into the primal value.

if sol.info.problem == 1 
    if strcmp(sol.objs,'inf')
        sol.pval = inf;
    else
        sol.pval = -inf;
    end
elseif sol.info.problem == 2
    if strcmp(sol.objs,'inf')
        sol.pval = -inf;
    else
        sol.pval = inf;
    end
end
end

function computedualsres(rel)

% Only PSD computations are working at the moment.

temp = 0;
if ~isempty(rel.F)
    temp = rel.F'*rel.lambda; 
end

if strcmp(rel.objs,'inf')
    temp = temp + rel.b;
else
    temp = temp - rel.b;
end

% Onto the support constraints.

rel.dres{1,1} = 'Equality constraints';

switch rel.rtyp
    case {'D','DD'}

        for j = 1:numel(rel.supcon)+1
            if isempty(rel.F{1})
                X{j} = dual(rel.ycons(j));
            else
                X{j} = dual(rel.ycons(2+j));
            end
        end

        for j = 1:numel(rel.supcon)+1
            temp = temp + rel.A{j}*X{j};
        end

        rel.dres{1,2} = -max(abs(temp));   % Return Linfty norm of equality constraing violations
        clear temp

        rel.dres{2,1} = 'Cone constraints';
        for j = 1:numel(rel.supcon)+1
            rel.dres{2,2} = [rel.dres{2,2};min(X{j})];
        end

    case 'SDD'

        if isempty(rel.mass) && isempty(rel.F{1})
            L = 0;
        elseif isempty(rel.mass) || isempty(rel.F{1})
            L = 1;
        else
            L = 2;
        end

        % Initialise SOC dual variables

        X{1} = zeros(3*nchoosek(n+d,d),1);

        for j = 1:numel(rel.supcon)
            dc = deg(rel.supcon(j));
            X{j+1} = zeros(3*nchoosek(n+floor(d-dc/2),floor(d-dc/2),1));
        end

        % Populate the SOC duals

        J = 1;
        for j = 1:numel(rel.supcon)+1
            temp2 = rel.A{j};
            for k = 1:numel(temp2)
                X{j}(1+3*(k-1):3*k,:) = dual(rel.ycons(L+J));
                temp = temp + temp2{k}*X{j}(1+3*(k-1):3*k,:);
                J = J + 1;
            end
            clear temp2
        end

        rel.dres{1,2} = -max(abs(temp));    % Return Linfty norm of equality constraing violations
        clear temp

        rel.dres{2,1} = 'Cone constraints';
        rel.dres{2,2} = [];
        for j = 1:numel(rel.supcon)+1
            temp = X{j};
            rel.dres{2,2} = [rel.dres{2,2}; X{j}(1)-sqrt(X{j}(2)^2+X{j}(3)^2)];
            for k = 2:numel(X{j})/3
                rel.dres{2,2} = min(rel.dres{2,2},X{j}((k-1)*3+1)-sqrt(X{j}((k-1)*3+2)^2+X{j}((k-1)*3+3)^2));
            end
        end

    case 'FKW'

    case 'PSD'


        for j = 1:numel(rel.A)
            temp2 = [];
            A = rel.A{j};
            for k = 1:numel(rel.A{j})
                temp2 = [temp2;sum(sum(rel.X{j}.*A{k}))];
            end
            clear A;

            temp = temp + temp2; clear temp2
        end

        rel.dres{1,2} = -max(abs(temp)); % Return Linfty norm of equality constraing violations
        clear temp

        rel.dres{2,2} = [];

        for j = 1:numel(rel.A)
            rel.dres{2,1} = 'Cone constraints';
            rel.dres{2,2} = [rel.dres{2,2};min(eig(rel.X{j}))];
        end
end

end

