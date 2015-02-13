function out = shift(p,y)

% Shift sequence y by p.

% REMEMBER BECAUSE OF HOW sdpvar WORKS, THE OUTPUT MUST NEVER BE A SINGLE
% CHARACTER (OTHERWISE WE ACCIDENTLY START REMAINING sdpvars). DO A TEST
% FOR THIS USING VAROUT.

% Juan Kuntz, 13/02/2015

for i = 1:numel(p.coef(1,:))
    temp(:,i) = grlext(p.nvar,p.coef(2,i),p.choose);
end

% If y is sdpvar object
    
if isa(y,'sdpvar')
    
    tab = ncktab(p.nvar+p.deg+sum(grlex(p.nvar,numel(y))));
    
    for i=1:numel(y)
        tempi = grlext(p.nvar,i,tab);
        test = 0;
        for j = 1:numel(p.coef(1,:))
            temp2 = igrlext(tempi+temp(:,j),tab);
            if temp2 <= numel(y) 
                if test == 0
                    test = 1;
                    out(i,1) =  p.coef(1,j)*y(temp2);
                else
                    out(i,1) =  p.coef(1,j)*y(temp2) + out(i);
                end
            end
        end
        clear tempi
    end
    return
end

% Else y is a seq object.

out = seq(p.nvar,p.deg+y.ord); % Initialise seq of appropiate dimension.
tab = ncktab(p.nvar+p.deg+y.ord);

ARRAY = y.coef; % This is needed because subsref for seq objects is not well developed (for example, it cannot handle y.coef(1,1) if y is seq).
ARRAY2 = [];

for i = 1:numel(ARRAY(1,:))
    tempi = grlext(p.nvar,i,tab);
    test = 0;
    for j = 1:numel(p.coef(1,:))
        temp2 = igrlext(tempi+temp(:,j),tab);
        I = bfind(ARRAY(2,:),temp2);
        if ~isempty(I)
            if test == 0
                test = 1;
                ARRAY2 = [ARRAY2,[p.coef(1,j)*ARRAY(1,I);i]];
            else
                ARRAY2(1,end) = ARRAY2(1,end) + p.coef(1,j)*ARRAY(1,I);
            end
        end
    end
    clear tempi
end

out.coef = ARRAY2;

end