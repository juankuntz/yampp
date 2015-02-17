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

% Actually solve
switch cp.reltype
    case 'D'
        
    case 'DD'
        
    case 'SDD'
        
    case 'FKW'
        
    case 'PSD'
        for i = 1:numel(cp.obj)
            temp = optimize(cp.ycons,cp.yobj(i),cp.ops); % Compute solution.
            
            % Store solution.
            
            cp.sol{end+1}.reltype = cp.reltype;
            cp.sol{end}.relorder = cp.relorder;
            cp.sol{end}.obj = cp.obj(i);
            cp.sol{end}.pval = double(cp.yobj(i));
            cp.sol{end}.ppoint = seq(cp.nvar,2*cp.relorder,double(cp.yvar));
            [cp.sol{end}.pres,~] = check(cp.ycons);
            cp.sol{end}.info = temp; clear temp;
            
            cp.sol{end}.dval = double(dual(cp.ycons(1)));
            
            for j = 1:numel(cp.seqeqcon)+1
                t(j,1) = dual(cp.ycons(j));
            end
            
            cp.sol{end}.dres = cp.F'*t-cp.b{i}; 
            clear t;
            
            for j = 1:numel(cp.supcon)+1
                X{j} = dual(cp.ycons(numel(cp.seqeqcon)+1+j));
            end
            
            for j = 1:numel(cp.supcon)+1
                temp = [];
                A = cp.A{j};
                for k = 1:numel(cp.A{j})
                    temp = [temp;sum(sum(X{j}.*A{k}))];
                end
                clear A;
            
                cp.sol{end}.dres = cp.sol{end}.dres + temp;
                
                % ADD THE RESIDUES OF THE PSD CONSTRAINTS, X{i}>=0 for all
                % i.
            end
        end
        
    case 'ALL'
        
end

end
