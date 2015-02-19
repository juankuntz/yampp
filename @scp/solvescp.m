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
        temp = optimize(cp.ycons,cp.yobj(i),cp.ops); % Compute solution.

        % Store solution.
        
        cp.sol{end+1}.reltype = cp.reltype;
        
        cp.sol{end}.ycons = cp.ycons;
        cp.sol{end}.relorder = d;
        cp.sol{end}.obj = cp.obj(i);
        cp.sol{end}.pval = value(cp.yobj(i));
        cp.sol{end}.ppoint = seq(n,2*d,value(cp.yvar));
        [cp.sol{end}.pres,~] = check(cp.ycons);
        cp.sol{end}.info = temp; clear temp;

        cp.sol{end}.dval = -value(dual(cp.ycons(1))); % Not sure why we need the minus sign here...

        
        t = dual(cp.ycons(1)); % Dual of the mass constraint.
        temp = dual(cp.ycons(2)); % Dual of the equality constraints
        t(2:numel(temp)+1,1) = temp;
        clear temp
        
        cp.sol{end}.dres = cp.F'*t; 
        clear t;

            
    switch cp.reltype
        case 'D' 
            
            % This is strange, but it seems we need to switch the sign here
            % depending on whether it's an LP, SOCP, or SDP...
            
            cp.sol{end}.dres = cp.sol{end}.dres + cp.b{i};
            
            % Onto the support constraints.
            
            for j = 1:numel(cp.supcon)+1
                X{j} = dual(cp.ycons(2+j));
            end
            
            for j = 1:numel(cp.supcon)+1
                cp.sol{end}.dres = cp.sol{end}.dres + cp.A{j}*X{j};
            end
            
            for j = 1:numel(cp.supcon)+1
                cp.sol{end}.dres = [cp.sol{end}.dres;X{j}];
            end
            
        case 'DD'
            
        case 'SDD'

        case 'FKW'

        case 'PSD'
            
            % This is strange, but it seems we need to switch the sign here
            % depending on whether it's an LP, SOCP, or SDP...
            
            cp.sol{end}.dres = cp.sol{end}.dres + cp.b{i};
            
                for j = 1:numel(cp.supcon)+1
                    X{j} = dual(cp.ycons(2+j));
                end

                for j = 1:numel(cp.supcon)+1
                    temp = [];
                    A = cp.A{j};
                    for k = 1:numel(cp.A{j})
                        temp = [temp;sum(sum(X{j}.*A{k}))];
                    end
                    clear A;

                    cp.sol{end}.dres = cp.sol{end}.dres + temp;
                end

                for j = 1:numel(cp.supcon)+1
                    cp.sol{end}.dres = [cp.sol{end}.dres;min(eig(X{j}))];
                end
        case 'ALL'
    end
        
end

end
