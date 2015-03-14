function out = getsols(cp,C)

% Find solutions that match the description stored in C. C must be an
% array, each element of which may specify the objective functions we are
% interested in, whether we are interested in both maxisimiations and
% minimisations, the relaxation degrees, and the type of relaxations. For
% example, C = {PSD,4:7,D,[x(1),x(2)],sup} tells getsols to return solutions
% of relexations types PSD or D, with objective function x(1) or x(2), of
% relaxation order 4,5,6, or 7, and that were maximisations. 

types = []; d = []; objf = []; objs = []; out = [];

% Extract and sort solution descriptions from C. 

for i = 1:numel(C)
    if isa(C{i},'double')
        
        d = [d;C{i}(:)];
        
    elseif isa(C{i},'pol')
        
        objf = [objf,C{i}(:)];
        
    elseif strcmpi(C{i},'inf') || strcmpi(C{i},'sup') 
        
        flg = 1;
        for j = 1:numel(objs)
            if strcmpi(objs{j},C{i})
                flg = 0;
            end
        end
        
        if flg
            objs{end+1} = C{i};
        end
        
    elseif strcmpi(C{i},'D') || strcmpi(C{i},'DD')  || strcmpi(C{i},'SDD')  || strcmpi(C{i},'FWK')  || strcmpi(C{i},'PSD') 
        
        flg = 1;
        for j = 1:numel(types)
            if strcmpi(types{j},C{i})
                flg = 0;
            end
        end
        
        if flg
            types{end+1} = C{i};
        end
    end
end

% Find solutions that match the descriptions and return them.

for i = 1:numel(cp.sol)
    test = 1;
    if ~isempty(d) && any(cp.sol{i}.relorder == d) == 0
        test = 0;
    elseif ~isempty(objf) && any(cp.sol{i}.obj == objf) == 0
        test = 0;
    elseif ~isempty(objs) 
        flg = 1;
        for j = 1:numel(objs)
            if strcmpi(objs{j},cp.sol{i}.minmax)
                flg = 0;
            end
        end
        if flg
            test = 0;
        end
    elseif ~isempty(types)
        flg = 1;
        for j = 1:numel(types)
            if strcmpi(types{j},cp.sol{i}.reltype)
                flg = 0;
            end
        end
        if flg
            test = 0;
        end
    end
    if test
        out{end+1} = cp.sol{i};
    end
end

end