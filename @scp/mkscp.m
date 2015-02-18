function mkscp(cp)

% Construct the appropiate moment problem.

% ADD WARNINGS/COMPLAINS IF HAVE ALREADY BEEN SOLVED IN THE PAST (SO THAT
% THE USER IS NOT UNDER THE IMPRESIONS THAT ALL SOLUTIONS CORRESPOND TO THE
% SAME DATA).

% Juan Kuntz, 13/02/2015, last edited 17/02/2015.

% Declare shorthands.

n = cp.nvar; d = cp.relorder;

% Clear any previous yalmip constraints.

cp.ycons = []; cp.A = []; cp.b = []; cp.F = []; cp.f = []; % LATER MAKE OPTION THAT SKIPS RE-WRITING b F and f if we are just switching type of constraints.

% Initialise yalmip variables.

cp.yvar = [];
cp.yvar = sdpvar(nchoosek(n+2*d,2*d),1); 

y = cp.yvar;

% Add mass constraints

if ~isempty(cp.mass)
    cp.ycons = [cp.ycons,y(1) - cp.mass == 0];
    cp.F = [1,zeros(1,nchoosek(n+2*d,2*d)-1)];   % Required later to compute dual residues.
    cp.f = 1;
end

% Add equality constraints.

for i = 1:numel(cp.seqeqcon)
    temp = coefficients(cp.seqeqcon(i))';
    cp.F = [cp.F; [temp,zeros(1,nchoosek(n+2*d,2*d)-numel(temp))]]; % Required later to compute dual residues.
    clear temp
    cp.f = [cp.f;0];
end

cp.ycons = [cp.ycons,cp.F(2:end,:)*y == 0];

% Add (linear) inequality constraints NOT IMPLEMENTED: CAREFUL WHEN
% IMPLEMENTING, IT AFFECTS SOME INDEXING IN SOLVESCP.M. FOR EXAMPLE IN
% switch cp.reltype
%        case 'D'
%            for j = 1:numel(cp.supcon)+1
%                X{j} = dual(cp.ycons(2+j));
%            end
% WE WOULD HAVE TO SWAP dual(cp.ycons(3+j)) for dual(cp.ycons(2+j));

% for i = 1:numel(cp.seqineqcon)
%     cp.ycons = [cp.ycons,cp.seqineqcon(i)*y >= 0];
% end

% Add objective

cp.yobj = [];

for i = 1:numel(cp.obj)
    cp.yobj = [cp.yobj,cp.obj(i)*y];
    
    % Store for use later retrieving the dual residuals.
    
    temp = coefficients(cp.obj(i));
    cp.b{i} = [temp;zeros(nchoosek(n+2*d,2*d)-numel(temp),1)];
    clear temp
end

% Add moment constraints.

switch cp.reltype
    case 'D'
        di(cp);
    case 'DD'
        
    case 'SDD'
        
    case 'FKW'
        
    case 'PSD'
        psd(cp);
end

end


function di(cp) % D constraints.

% Declare shorthands.

n = cp.nvar; d = cp.relorder;
y = cp.yvar;

tab = ncktab(n+2*d);

% General moment constraint

T = zeros(nchoosek(n+d,d),nchoosek(n+2*d,2*d));

for j = 1:nchoosek(n+d,d)
    T(j,igrlext(2*grlext(n,j,tab),tab)) = 1;
    
end
       
cp.ycons = [cp.ycons,T*y >= 0];
cp.A{1} = T'; 

% Localising constraints

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    l(1) = nchoosek(n+floor(d-cp.supcon(i).deg/2),floor(d-cp.supcon(i).deg/2));
    l(2) = nchoosek(n+2*d-cp.supcon(i).deg,2*d-cp.supcon(i).deg);
    
    [temp,Tq] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    
    cp.ycons = [cp.ycons,T(1:l(1),1:l(2))*temp >= 0];
    cp.A{end+1} = (T(1:l(1),1:l(2))*Tq)';
    
    clear temp Tq l
end

end

function psd(cp) % PSD constraints.

% Declare shorthands.

n = cp.nvar; d = cp.relorder;
y = cp.yvar;



% Moment matrix constraint.

B = hankelbasis(n,d);

% We need this when we recover the dual residues later.

for i = 1:numel(B)
    temp{i} = - B{i};
end
cp.A{1} = temp;    
clear temp

Ay = zeros(size(B{1}));
for i = 1:nchoosek(n+2*d,2*d)
    Ay = Ay + B{i}*y(i);
end

cp.ycons = [cp.ycons,Ay >= 0];

clear B; clear Ay;

% Localising matrices constraints.

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    % Construct localising matrix.
    
    B = hankelbasis(n,floor(d-cp.supcon(i).deg/2));
    [temp,T] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    
    % We need the follwing when we recover the dual residues later.
    
    for j = 1:numel(cp.yvar) 
        temp2{j} = zeros(size(B{1}));
        for k = 1:numel(B)
            Tt = T';
            temp2{j} = temp2{j} - Tt(j,k)*B{k};
        end
    end
    
    cp.A{end+1} = temp2;    
    clear temp2

    % Back to constructing the locasing matrix.
    
    Ay = zeros(size(B{1}));    
    
    for j = 1:nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2))
        Ay = Ay + B{j}*temp(j);
    end
    
    cp.ycons = [cp.ycons,Ay >= 0]; % Add localising matrix constraint.
    
    clear temp Ay
end

end

