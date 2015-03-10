function out = disp(p)

% Displays the polynomial.

% Juan Kuntz, 02/02/2015

[n,m] = size(p);
out{n,m} = [];
for i = 1:n
    for j = 1:m
        out{i,j} = computestring(p(i,j)); % Do it entrywise
    end
end
disp(out)

end

function out = computestring(p)

% Construct the output string for the scalar polynomial p.

if isempty(p.coef) % Constant polynomial
    out = '0';
    return
elseif numel(p.coef(1,:)) == 1 && p.coef(2,1) == 1
    out = num2str(p.coef(1,1));
    return
end

% Non-trivial polynomial

FIRST = 1; out = ''; n = sum(p.var.ncomp);

TAB = ncktab(p.nvar+p.deg); % Avoids repeated calls to ncktab.
for i = 1:numel(p.coef(1,:))
    mono = grlext(n,p.coef(2,i),TAB);
    if p.coef(1,i)<0
        out =  strcat(out,'-');
    end
    if FIRST
        FIRST = 0;
        if abs(p.coef(1,i)) ~= 1 || sum(mono) == 0
            out =  strcat(num2str(p.coef(1,i)));
        end
    else
        if p.coef(1,i) >= 0
            out =  strcat(out,'+');
        end
        if abs(p.coef(1,i)) ~= 1 || sum(mono) == 0
            out =  strcat(out,num2str(abs(p.coef(1,i))));
        end
    end
    temp = 1; l = 1;
    for k = 1:n
       if mono(k) == 1
        out = strcat(out,p.var.symb(temp),'(',num2str(l),')');
       elseif mono(k) ~= 0
           out = strcat(out,p.var.symb(temp),'(',num2str(l),')^',num2str(mono(k)));
       end
       if k == sum(p.var.ncomp(1:temp)) 
           temp = temp + 1; l = 1;
       else 
           l = l + 1;
       end
    end
    clear mono
end

end

