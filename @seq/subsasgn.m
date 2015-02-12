function y = subsasgn(y,S,c)

% Check for no repeated monomials and also that length(c) == l also that
% all monomial vectors are the same length (that of the dimension of the
% thingy)

% Explain the three types of indexing...

% Possible improvement in speed if argument is given sorted?

% Juan Kuntz, 12/02/2015

% Allow access to properties.

if strcmp(S.type,'.') 
    switch S.subs
        case 'coef'
            y.coef = c;
        case 'dim'
            y.dim = c;
        case 'ord'
            y.ord = c;
        case 'choose'
            y.choose = c;
    end
    return
end

% Direct access to sequence element mon.

if y.dim > 1 && numel(S.subs) == y.dim 
    for i = 1:numel(S.subs) % Extract indexes.
        mon(i) = S.subs{i};
    end
    if isempty(y.ord) || y.ord < sum(mon) % Update order and choose table need be.
        y.ord = sum(mon);
        y.choose = ncktab(sum(mon)+y.dim);
    end
    rank = igrlext(mon,y.choose);
    if isempty(y.coef)
        y.coef = [c;rank];
    else
        I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
        if isempty(I)
            I = bfindint(y.coef(2,:),rank);
            y.coef = [y.coef(:,1:I),[c;rank],y.coef(:,I+1:end)];
        else
            y.coef(1,I) = c;
        end
    end
    return
end

if numel(c) ~= numel(S.subs{1}(1,:))
    disp('Error: The number of assignemnts must be the same as that of assignees,');
end

mon = S.subs{1};

% Coefficients specified by a matrix whose columns are exponents.

if isdouble(mon)
    % MISSING: Check that there are no repeated indexes. ALSO THAT ALL
    % EXPONENTS ARE VALID
    for i = 1:numel(c) % SEARCHING HERE IS SUBOPTIMAL
        if isempty(y.ord) || y.ord < sum(mon(:,i)) % Update order and choose table need be.
            y.ord = sum(mon(:,i));
            y.choose = ncktab(sum(mon(:,i))+y.dim);
        end
        rank = igrlext(mon(:,i),y.choose);
        if isempty(y.coef)
            y.coef = [c(i);rank];
        else
            I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
            if isempty(I)
                I = bfindint(y.coef(2,:),rank);
                y.coef = [y.coef(:,1:I),[c(i);rank],y.coef(:,I+1:end)];
            else
                y.coef(1,I) = c(i);
            end
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

for i = 1:numel(c) % SEARCHING HERE IS SUBOPTIMAL
    if isempty(y.ord) || y.ord < mon(i).deg % Update order and choose table need be.
        y.ord = mon(i).deg;
        y.choose = ncktab(mon(i).deg+y.dim);
    end
    rank = mon(i).coef(2);
    if isempty(y.coef)
        y.coef = [c(i);rank];
    else
        I = bfind(y.coef(2,:),rank); % SEARCHING HERE IS SUBOPTIMAL
        if isempty(I)
            I = bfindint(y.coef(2,:),rank);
            y.coef = [y.coef(:,1:I),[c(i);rank],y.coef(:,I+1:end)];
        else
            y.coef(1,I) = c(i);
        end
    end
end
    

end