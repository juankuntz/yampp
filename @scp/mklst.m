function mklst(cp)

% Strings out the basic information regarding each relaxation to be solved.

% The field names are different because of historical reasons, they should
% be updated.

% Juan Kuntz, 14/03/2015, last edited 12/05/2015.

objs = cp.obj{1}; objf = cp.obj{2};

d = cp.rord; typ = cp.rtyp;

out = [];
for i = 1:numel(typ)
    for j = 1:numel(d)
        for k = 1:numel(objf)
            
            out{end+1}.rord =  d(j);
            out{end}.mult = cp.mult^cp.multpow(d(j));
            out{end}.objf = objf(k);
            out{end}.objs = objs(k,:);
            
            out{end}.FW = [];
            if ischar(typ{i}) == 1
                out{end}.rtyp = typ{i};
            else
                temp = typ{i};
                out{end}.rtyp = temp{1};
                out{end}.FW = temp{2};
            end
            
            % Fish out solver options.

            switch typ{i}
                case {'d','D'}
                    out{end}.ops = cp.ops{1,end};
                case {'dd','DD'}
                    out{end}.ops = cp.ops{2,end};
                case {'sdd','SDD'}
                    out{end}.ops = cp.ops{3,end};
                case {'fwk','FWK'}
                    out{end}.ops = cp.ops{4,end};
                case {'psd','PSD'}
                    out{end}.ops = cp.ops{5,end};
                case {'nn','NN'}
                    out{end}.ops = cp.ops{6,end};
            end

        end
    end
end
cp.rlst = out;
end
