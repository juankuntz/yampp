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
    cp.ycons = [cp.ycons,y(1) == cp.mass];
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
        dd(cp);
    case 'SDD'
        sdd(cp);
    case 'FWK'
        fwk(cp);
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
cp.A{1} = -T'; 

% Localising constraints

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    l(1) = nchoosek(n+floor(d-cp.supcon(i).deg/2),floor(d-cp.supcon(i).deg/2));
    l(2) = nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2));
    
    [temp2,Tq2] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l(2)); 
    Tq = Tq2(1:l(2),:);
    
    clear temp2 Tq2
    
    cp.ycons = [cp.ycons,T(1:l(1),1:l(2))*temp >= 0];
    cp.A{end+1} = -(T(1:l(1),1:l(2))*Tq)';
    
    clear temp Tq l
end

end

function dd(cp) % DD constraints.

% Declare shorthands.

n = cp.nvar; d = cp.relorder;
y = cp.yvar;

tab = ncktab(n+2*d);

% General moment constraint

l = nchoosek(n+2*d,2*d);
T = zeros(nchoosek(n+d,d),l);

mon = zeros(cp.nvar,l);

% Diagonal entries of the moment matrix must be greater or equal than zero.

for i = 1:nchoosek(n+d,d) 
    mon(:,i) = grlext(n,i,tab);
    T(i,igrlext(2*mon(:,i),tab)) = 1; %y_2a >= 0 
end

% Now the off-diagonal constaints.

f = T*(1:l)';
for i = 1:nchoosek(n+d,d)
    for j = 1:i-1
        T = [T;zeros(2,l)];
        k = igrlext(mon(:,i)+mon(:,j),tab);
        T(end-1,f(i)) = 1; T(end-1,f(j)) = 1; T(end-1,k) = -2; % y_2a+y_2b - 2y_(a+b)>=0
        T(end,f(i)) = 1; T(end,f(j)) = 1; T(end,k) = 2; % y_2a+y_2b - 2y_(a+b)>=0
    end
end

clear f l mon k

cp.ycons = [cp.ycons,T*y >= 0];
cp.A{1} = -T'; 

% Localising constraints

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    l(1) = nchoosek(n+floor(d-cp.supcon(i).deg/2),floor(d-cp.supcon(i).deg/2));
    l(2) = nchoosek(n+d,d);
    l(3) = 0;
    for j = 1:l(1)
        l(3) = l(3) + 2*(j-1);
    end
    l(4) = nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2));
    
    
    [temp2,Tq2] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l(4));
    Tq = Tq2(1:l(4),:);
    
    clear temp2 Tq2
    
    C = [T(1:l(1),1:l(4));              % Diagonal constaints
        T(l(2)+1:l(2)+l(3),1:l(4))];  % Off diagonal constaints.
    
    cp.ycons = [cp.ycons,C*temp >= 0]; 
    
    cp.A{end+1} = -(C*Tq)';
    
    clear temp Tq l C
end


end

function sdd(cp) % SDD constraints.

% Declare shorthands.

n = cp.nvar; d = cp.relorder;
y = cp.yvar; l2d = nchoosek(n+2*d,2*d);

tab = ncktab(n+2*d);

% General moment constraint

mon = zeros(cp.nvar,l2d);
r = zeros(cp.nvar,1);

for i = 1:nchoosek(n+d,d) 
    mon(:,i) = grlext(n,i,tab);
    r(i) = igrlext(2*mon(:,i),tab);
end

l = 1;
for i = 1:nchoosek(n+d,d)
    for j = 1:i-1
        
        k = igrlext(mon(:,i)+mon(:,j),tab);
        cp.ycons = [cp.ycons,cone([2*y(k);y(r(i))-y(r(j))],y(r(i))+y(r(j)))];
        
        T = zeros(3,l2d);
        T(1,r(i)) = 1; T(1,r(j)) = 1;
        T(2,k) = 2;
        T(3,r(i)) = 1; T(3,r(j)) = -1;
        
        Tt{l} = -T';
        
        l = l + 1;
    end
end

cp.A{1} = Tt;  

clear T Tt

% Localising matrices constraints.

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    l2dn = nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2));
    
    [temp2,Tq2] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l2dn);
    Tq = Tq2(1:l2dn,:);
    
    l = 1;
    for j = 1:nchoosek(n+floor(d-cp.supcon(i).deg/2),floor(d-cp.supcon(i).deg/2))
        for k = 1:j-1

            I = igrlext(mon(:,j)+mon(:,k),tab);
            cp.ycons = [cp.ycons,cone([2*temp(I);temp(r(j))-temp(r(k))],temp(r(j))+temp(r(k)))];

            T = zeros(3,l2dn);
            T(1,r(j)) = 1; T(1,r(k)) = 1;
            T(2,I) = 2;
            T(3,r(j)) = 1; T(3,r(k)) = -1;
            
            Tt{l} = -Tq'*T';
            
            l = l + 1;
        end
    end

    cp.A{end+1} = Tt;  
    clear T Tt
end

end

function fwk(cp)

% Declare shorthands.

n = cp.nvar; d = cp.relorder;
y = cp.yvar; FW = cp.FW;

% Moment matrix constraint.

B = hankelbasis(n,d);

M = zeros(size(B{1}));
for i = 1:nchoosek(n+2*d,2*d)
    M = M - (-B{i})*y(i);
end

l = nchoosek(n+d,d);
I = nchoosek(1:l,FW);

for j = 1:numel(I(:,1))
    cp.ycons = [cp.ycons,M(I(j,:),I(j,:)) >= 0];
end

clear temp Ay l I

%% THIS SHOULD WORK (AND REDUCE THE SPEED ISSUES) BUT YALMIP RUNS OUT OF MEMORY.

% Ays = []; X = []; Y = [];
% for i = 1:numel(I(:,1))
%     Ays = [Ays;vec(M(I(i,:),I(i,:)))];
%     X = [X;(i-1)*FW*ones(FW^2,1)+vec(repmat((1:FW)',[1,FW]))];
%     Y = [Y;(i-1)*FW*ones(FW^2,1)+vec(repmat(1:FW,[FW,1]))];
% end
% 
% cp.ycons = [cp.ycons,sparse(X,Y,Ays,numel(I(:,1))*5,numel(I(:,1))*5) >= 0];

clear l I Ay X Y

% Localising matrices constraints.

for i = 1:numel(cp.supcon)
    
    if 2*d < cp.supcon(i).deg
        disp(['Error: the degree of the ',num2str(i),'"s polynomial defining the support is too big']); % fix this.
        return
    end
    
    % Construct localising matrix.
    
    B = hankelbasis(n,floor(d-cp.supcon(i).deg/2));
    [temp,T] = shift(cp.supcon(i),y);   % Shift the sequence by the support polynomial.
    
    % Back to constructing the localising matrix.
    
    Ay = zeros(size(B{1}));    
    
    for j = 1:nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2))
        Ay = Ay - (-B{j})*temp(j);
    end

    l = nchoosek(n+floor(d-cp.supcon(i).deg/2),floor(d-cp.supcon(i).deg/2));
    I = nchoosek(1:l,FW);

    if l < FW
        cp.ycons = [cp.ycons,Ay >= 0];
    else
        for j = 1:numel(I(:,1))
            cp.ycons = [cp.ycons,Ay(I(j,:),I(j,:)) >= 0];
        end
    end
    clear temp Ay l I
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
    Ay = Ay - (-B{i})*y(i);
end

cp.ycons = [cp.ycons,Ay >= 0]; % C - A^T(y) >= 0

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

    % Back to constructing the localising matrix.
    
    Ay = zeros(size(B{1}));    
    
    for j = 1:nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2))
        Ay = Ay - (-B{j})*temp(j);
    end
    
    cp.ycons = [cp.ycons,Ay >= 0]; % Add localising matrix constraint.
    
    clear temp Ay
end

end

