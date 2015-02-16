function [pnotq,qnotp] = varcomp(p,q)

% Finds which variable symbols are part of p but not part of q, and stores
% theses in pnotq (along with their number of components). Ditto for qnotp.

% Juan Kuntz, 06/02/2015

pnotq.symb = []; pnotq.ncomp = [];
qnotp.symb = []; qnotp.ncomp = [];

both = [];

for i = 1:numel(p.var.symb)
    for j = 1:numel(q.var.symb)
        if p.var.symb(i) == q.var.symb(j)
            both = [both,p.var.symb(i)];
            break
        end
    end
end

for i = 1:numel(p.var.symb)
    test = 0;
    for j = 1:numel(both)
        if p.var.symb(i) == both(j)
            test = 1;
        end
    end
    if test == 0
        pnotq.symb = [pnotq.symb,p.var.symb(i)];
        pnotq.ncomp = [pnotq.ncomp,p.var.ncomp(i)];
    end
end


for i = 1:length(q.var.symb)
    test = 0;
    for j = 1:numel(both)
        if q.var.symb(i) == both(j)
            test = 1;
        end
    end
    if test == 0
        qnotp.symb = [qnotp.symb,q.var.symb(i)];
        qnotp.ncomp = [qnotp.ncomp,q.var.ncomp(i)];
    end
end

end


