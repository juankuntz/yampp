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
    if ~isempty(cp.eqcon{1})
        sol.dval = -sol.f'*value(dual(sol.ycons(1))); % Minus sign here because yalmip fits this to a maximisation problem, see Lofberg's paper on Dualize.
    end
elseif strcmp(sol.objs,'sup')
    temp = optimize(sol.ycons,-sol.yobj,sol.ops);
    if ~isempty(cp.eqcon{1})
        sol.dval = sol.f'*value(dual(sol.ycons(1))); 
    end
end

% Store solution.

sol.pval = value(sol.yobj);
sol.ppoint = seq(n,2*d,value(sol.yvar));

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
