function out = subsref(y,S)

if strcmp(S.type,'.') % Allow access to properties.
    switch S.subs
        case 'coef'
            out = y.coef;
        case 'dim'
            out = y.dim;
        case 'ord'
            out = y.ord;
        case 'choose'
            out = y.choose;
    end
    return
end

% Direct access to sequence element mon.

if y.dim > 1 && numel(S.subs) == y.dim 
    for i = 1:numel(S.subs) % Extract indexes.
        mon(i) = S.subs{i};
    end
    rank = igrlext(mon,y.choose);
    I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
    if isempty(I)
        out = 0;
    else
        out = y.coef(1,I);
    end
    return
end

mon = S.subs{1};
out = [];

% Coefficients specified by a matrix whose columns are exponents.

if isdouble(mon)
    % MISSING: Check that there are no repeated indexes. ALSO THAT ALL
    % EXPONENTS ARE VALID
    for i = 1:numel(mon(1,:)) % SEARCHING HERE IS SUBOPTIMAL
        rank = igrlext(mon(:,i),y.choose);

        I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
        if isempty(I)
            out = [out,0];
        else
            out = [out,y.coef(1,I)];
        end
        
    end
    return
end

% Final option, coefficients specified by a vector of monomials: First 
% check that all indexes are actually valid (that is monomials).

test = 1;
for i = 1:numel(mon)
    if ~ismon(mon(i))
        test = 0;
    end
end

if test == 0
    disp('The only polnoymials tha can be used  monomials can be used as indexes.');
    return
end

% Actually do the assignments.

for i = 1:numel(mon) % SEARCHING HERE IS SUBOPTIMAL
    rank = mon(i).coef(2);
    I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
    if isempty(I)
        I = bfindint(y.coef(2,:),rank);
        out = [out,0];
    else
        out = [out,y.coef(1,I)];
    end
end

end