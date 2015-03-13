function mklst(cp)

% Strings out the basic information regarding each relaxation to be solved.

% The field names are different because of historical reasons, they should
% be updated.

% Juan Kuntz, 13/03/2015.

objs = cp.obj{1}; objf = cp.obj{2};

d = cp.rord; typ = cp.rtyp;

out = [];
for i = 1:numel(typ)
    for j = 1:numel(d)
        for k = 1:numel(objf)
            
            out{end+1}.nvar =  cp.nvar;
            out{end}.mass = cp.mass;
            out{end}.seqeqcon = cp.eqcon{1};
            out{end}.ineqcon = cp.ineqcon{1};
            out{end}.supcon = cp.supineq;
            out{end}.ops = cp.ops;
            %out. = cp.supeqcon;
            
            
            out{end}.relorder =  d(j);
            out{end}.obj = objf(k);
            out{end}.minmax = objs(k,:);
            
            out{end}.FW = [];
            if ischar(typ{i}) == 1
                out{end}.reltype = typ{i};
            else
                temp = typ{i};
                out{end}.reltype = temp{1};
                out{end}.FW = temp{2};
            end

        end
    end
end
    
cp.rlst = out;
end
