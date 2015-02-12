function out = summi(n,d)

% [j;k] is a column of out{i} if and only if 
% grlex(n,j) + grlex(n,k) = grlex(n,i).
% 
% The function finds the above pairs for a deg(grlex(i)) <= 2d and
% deg(grlex(j)), deg(grlex(k)) <= d.
%
% Juan Kuntz, 12/02/2015.

tab = ncktab(n+2*d);

out{tab(n+2*d+1,2*d+1)} = [];

monnd = zeros(n,tab(n+d+1,d+1));

for i = 1:tab(n+d+1,d+1)
    monnd(:,i) = grlext(n,i,tab);
end

for i = 1:tab(n+d+1,d+1)
    for j = 1:i
        temp = igrlext(monnd(:,i)+monnd(:,j),tab);
        out{temp} = [out{temp},[i;j]];
        if i ~= j
            out{temp} = [out{temp},[j;i]];
        end
        clear temp
    end
end

end

