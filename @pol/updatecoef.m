function obj = updatecoef(obj,nvarsold,nsymbnew,oldsymb,oldchoose)

% If there is a change in the order and/or number variables of the polynomial obj, this function
% re-orders the coefficients of obj in the new grlex order. This function
% is only called by the set.var method.

% Juan Kuntz, 06/02/2015

if isempty(obj.coef) % If its the zero polynomial we are done, there is no re-ordering of coefficients needed (there are no coefficients to reoder!).
    return
end

% Otherwise, re-order the existing coefficients by the grlex order in
% the new variables (paying attention to the new
% alphabetical ordering of the symbols).

for i = 1:numel(obj.coef(1,:))
    oldmon = grlext(nvarsold,obj.coef(2,i),oldchoose); % Find the multiindex (in nvarsold number of variables) of the monomial corresponding to the coefficient obj.coef(2,i).

    % Now re-write this multiindex in the new variables
    % (using the correcting new alphabetical ordering).

    l = 0; u = 0;
    m = 1;
    for j = 1:nsymbnew
        test = 0;
        for k = 1:obj.var.ncomp(j)
            l = l + 1;
            if m > numel(oldsymb) || ~(oldsymb(m) == obj.var.symb(j)) % If we have not yet run out of variables, check whether the current variable is one of the old ones. If it is not, the multiindex entry is simply zero.
                newmon(l,1) = 0;
            else % If it is, then update according.
                u = u + 1;
                newmon(l,1) = oldmon(u);
                test = 1;
            end
        end
        if test == 1
                m = m + 1;
        end
    end
    newrank = igrlext(newmon,obj.choose);
    tempcoef(:,i) = [obj.coef(1,i);newrank];

    clear newmon newrank
end

obj.coef = tempcoef;
end

