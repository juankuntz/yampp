function solvescp(cp)

%
% Juan Kuntz, 16/02/2015, last edited 13/03/2015.

% Extract dual point and save in sol.
% 
% Compute, by hand (as in http://users.isy.liu.se/en/rt/johanl/2009_OMS_DUALIZE.pdf), the primal and dual residues, use dual() to extract dual variables.
% 
% Compare these with the residues returned by check().

% Missing error/warning messages.

setvars(cp);    % Decide the dimension of the underlying space.
mklst(cp);      % Make list of all relaxations to be solved.

% Build scp.

rels{numel(cp.rlst)} = []; sol{numel(cp.rlst)} = [];

if cp.parallelise == 1
    parfor i = 1:numel(cp.rlst)
        rels{i} = mkscp(cp.rlst{i});
        sol{i} = solve(rels{i});
    end
else
    for i = 1:numel(cp.rlst)
        rels{i} = mkscp(cp.rlst{i});
        sol{i} = solve(rels{i});
    end
end

for i = 1:numel(cp.rlst)
   cp.sol{i} = sol{i}; 
end

end


function sol = solve(cp)

% Shorthands

n = cp.nvar; d = cp.relorder;

%


if strcmp(cp.minmax,'inf')
    temp = optimize(cp.ycons,cp.yobj,cp.ops); % Compute solution.
    cp.sol.minmax = 'inf';

    cp.sol.dval = -value(dual(cp.ycons(1))); % Not sure why we need the minus sign here...
elseif strcmp(cp.minmax,'sup')
    temp = optimize(cp.ycons,-cp.yobj,cp.ops);
    cp.sol.minmax = 'sup';

    cp.sol.dval = value(dual(cp.ycons(1))); 
end

    % Store solution.

    cp.sol.reltype = cp.reltype;
    cp.sol.FW = cp.FW;
    cp.sol.obj = cp.obj;

    cp.sol.ycons = cp.ycons;
    cp.sol.relorder = d;

    cp.sol.pval = value(cp.yobj);
    cp.sol.ppoint = seq(n,2*d,value(cp.yvar));


    % Store primal residues in a formated table.

    [temp2,~] = check(cp.ycons);

    cp.sol.pres{1,2} = [];

    l = 0;
    if ~isempty(cp.mass)
        cp.sol.pres{1,1} = 'Mass constraint.';
        cp.sol.pres{1,2} = temp2(1);
        l = l + 1;
    end

    if ~isempty(cp.seqeqcon)
        if l == 0
            cp.sol.pres{1,1} = 'Equality constraints.';
            cp.sol.pres{1,2} = temp2(1);
        else
            cp.sol.pres{2,1} = 'Equality constraints.';
            cp.sol.pres{2,2} = temp2(2);
        end
        l = l + 1;
    end

    if l == 0
        cp.sol.pres{1,1} = 'Cone constraints.';
        cp.sol.pres{1,2} = temp2(1);
    else
        cp.sol.pres{end+1,1} = 'Cone constraints.';
        cp.sol.pres{end,2} = temp2(l+1);
    end
    l = l + 1;

   for j = 1:numel(cp.supcon)
        cp.sol.pres{end,2} = [cp.sol.pres{end,2};temp2(l+j)];
   end

    % Store yalmip output info.

    cp.sol.info = temp; clear temp temp2;

    % If problem is unbounded or infeasible store the appropiate +-inf
    % into the primal value.

    if cp.sol.info.problem == 1 
        if strcmp(cp.minmax,'inf')
            cp.sol.pval = inf;
        else
            cp.sol.pval = -inf;
        end
    elseif cp.sol.info.problem == 2
        if strcmp(cp.minmax,'inf')
            cp.sol.pval = -inf;
        else
            cp.sol.pval = inf;
        end
    end

    % If dual residues not requested, exit (saves considerable time).

    if cp.dualres == 0
        sol = cp.sol;
        return
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


    if strcmp(cp.minmax,'inf')
        temp = temp + cp.b;
    else
        temp = temp - cp.b;
    end

    % Onto the support constraints.

    cp.sol.dres{1,1} = 'Equality constraints';

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

            cp.sol.dres{1,2} = -max(abs(temp));   % Return Linfty norm of equality constraing violations
            clear temp

            cp.sol.dres{2,1} = 'Cone constraints';
            for j = 1:numel(cp.supcon)+1
                cp.sol.dres{2,2} = [cp.sol.dres{2,2};min(X{j})];
            end

        case 'SDD'

            if isempty(cp.mass) && isempty(cp.seqeqcon)
                L = 0;
            elseif isempty(cp.mass) || isempty(cp.seqeqcon)
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
                temp2 = cp.A{j};
                for k = 1:numel(temp2)
                    X{j}(1+3*(k-1):3*k,:) = dual(cp.ycons(L+J));
                    temp = temp + temp2{k}*X{j}(1+3*(k-1):3*k,:);
                    J = J + 1;
                end
                clear temp2
            end

            cp.sol.dres{1,2} = -max(abs(temp));    % Return Linfty norm of equality constraing violations
            clear temp

            cp.sol.dres{2,1} = 'Cone constraints';
            cp.sol.dres{2,2} = [];
            for j = 1:numel(cp.supcon)+1
                temp = X{j};
                cp.sol.dres{2,2} = [cp.sol.dres{2,2}; X{j}(1)-sqrt(X{j}(2)^2+X{j}(3)^2)];
                for k = 2:numel(X{j})/3
                    cp.sol.dres{2,2} = min(cp.sol.dres{2,2},X{j}((k-1)*3+1)-sqrt(X{j}((k-1)*3+2)^2+X{j}((k-1)*3+3)^2));
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

            cp.sol.dres{1,2} = -max(abs(temp)); % Return Linfty norm of equality constraing violations
            clear temp

            cp.sol.dres{2,2} = [];

            for j = 1:numel(cp.supcon)+1
                cp.sol.dres{2,1} = 'Cone constraints';
                cp.sol.dres{2,2} = [cp.sol.dres{2,2};min(eig(X{j}))];
            end
    end

    sol = cp.sol;
end
