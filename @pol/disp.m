function disp(p)

% Displays the polynomial.

% Juan Kuntz, 02/02/2015, last edited 15/03/2015.

[n,m] = size(p);
out{n,m} = [];
for i = 1:n
    for j = 1:m
        out{i,j} = computestring(p(i,j)); % Do it entrywise
    end
end

disp(out);

end


