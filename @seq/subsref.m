function out = subsref(y,S)

m = S.subs;

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

end