function rel = mkrel(cp,rspec)

% Construct the appropiate moment problem.

% Juan Kuntz, 13/02/2015, last edited 14/03/2015.

% Save data contained in rspec into rel and clear rspec.

rel = rspec;
clear rspec;

% Declare shorthands.

n = cp.nvar;
d = rel.relorder;

% THIS NEEDS FIXING IN THE FUTURE

rel.ops = cp.ops;
rel.dualres = 0;

% Initialise sdpvars and constraints

rel.yvar = sdpvar(nchoosek(n+2*d,2*d),1); 
y = rel.yvar;
rel.ycons = [];

% Add objective in terms of the sdpvars and store b vector in case we want 
% to compute the dual residues later.

rel.yobj = rel.obj*y;

temp = coefficients(rel.obj);
rel.b = [temp;zeros(nchoosek(n+2*d,2*d)-numel(temp),1)];
clear temp

% Add equality constraints.

rel.F = []; rel.f = [];
for i = 1:numel(cp.eqcon{2})
    if deg(cp.eqcon{1}(i)) <= 2*d
        temp = coefficients(cp.eqcon{1}(i))';
        rel.F = [rel.F; [temp,zeros(1,nchoosek(n+2*d,2*d)-numel(temp))]]; % Required later to compute dual residues.
        rel.f = [rel.f;cp.eqcon{2}(i)];
        clear temp  
    else
        warning(['At least one equality constraint has been ignored because the degree of that constraint`s polynomial is too large of the relaxations of order ',num2str(d),'.']);
    end
end

rel.ycons = [rel.ycons,rel.F*y == rel.f];

% Add (linear) inequality constraints NOT IMPLEMENTED: CAREFUL WHEN
% IMPLEMENTING, IT AFFECTS SOME INDEXING IN SOLVESCP.M. FOR EXAMPLE IN
% switch cp.reltype
%        case 'D'
%            for j = 1:numel(cp.supineq)+1
%                X{j} = dual(cp.ycons(2+j));
%            end
% WE WOULD HAVE TO SWAP dual(cp.ycons(3+j)) for dual(cp.ycons(2+j));

% for i = 1:numel(cp.seqineqcon)
%     cp.ycons = [cp.ycons,cp.seqineqcon(i)*y >= 0];
% end

% Add moment constraints.

switch rel.reltype
    case 'D'
        rel = di(cp,rel);
    case 'DD'
        rel = dd(cp,rel);
    case 'SDD'
        rel = sdd(cp,rel);
    case 'FWK'
        rel = fwk(cp,rel);
    case 'PSD'
        rel = psd(cp,rel);
end

end

function rel = di(cp,rel) % D constraints.

% Declare shorthands.

n = cp.nvar; 
d = rel.relorder;
y = rel.yvar;

ld = nchoosek(n+d,d);
l2d = nchoosek(n+2*d,2*d);
tab = ncktab(n+2*d);

% General moment constraint: Diagonal entries of the moment matrix must be 
% greater or equal than zero.

T = zeros(ld,l2d);

for j = 1:ld
    T(j,igrlext(2*grlext(n,j,tab),tab)) = 1;    
end
       
rel.ycons = [rel.ycons,T*y >= 0];
rel.A{1} = -T'; 

% Localising constraints

for i = 1:numel(cp.supineq)
    
    % Declare shorthands.
    
    dc = deg(cp.supineq(i));
    ldc = nchoosek(n+floor(d-dc/2),floor(d-dc/2));
    l2dc = nchoosek(n+2*floor(d-dc/2),2*floor(d-dc/2));
    
    [temp2,Tq2] = shift(cp.supineq(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l2dc); 
    Tq = Tq2(1:l2dc,:);
    
    clear temp2 Tq2
    
    rel.ycons = [rel.ycons,T(1:ldc,1:l2dc)*temp >= 0];
    rel.A{end+1} = -(T(1:ldc,1:l2dc)*Tq)';
    
    clear temp Tq dc ldc l2dc
end

end

function rel = dd(cp,rel) % DD constraints.

% Declare shorthands.

n = cp.nvar; 
d = rel.relorder;
y = rel.yvar;

ld = nchoosek(n+d,d);
l2d = nchoosek(n+2*d,2*d);
tab = ncktab(n+2*d);

% General moment constraint.

T = zeros(ld^2,l2d);
mon = zeros(n,ld);

% Diagonal entries of the moment matrix must be greater or equal than zero.

for i = 1:ld
    mon(:,i) = grlext(n,i,tab);
    T(i,igrlext(2*mon(:,i),tab)) = 1; %y_2a >= 0 
end

% Now the off-diagonal constaints.

u = ld;
f = T(1:nchoosek(n+d,d),:)*(1:l2d)';
for i = 1:nchoosek(n+d,d)
    for j = 1:i-1
        u = u + 2;
        k = igrlext(mon(:,i)+mon(:,j),tab);
        T(u-1,f(i)) = 1; T(u-1,f(j)) = 1; T(u-1,k) = -2; % y_2a+y_2b - 2y_(a+b)>=0
        T(u,f(i)) = 1; T(u,f(j)) = 1; T(u,k) = 2; % y_2a+y_2b - 2y_(a+b)>=0
    end
end

clear f mon k C

rel.ycons = [rel.ycons,T*y >= 0];
rel.A{1} = -T'; 

% Localising constraints

for i = 1:numel(cp.supineq)
    
    % Declare shorthands.
    
    dc = deg(cp.supineq(i));
    ldc = nchoosek(n+floor(d-dc/2),floor(d-dc/2)); 
    l2dc = nchoosek(n+2*floor(d-dc/2),2*floor(d-dc/2));
    
    k = 0;
    for j = 1:ldc
        k = k + 2*(j-1);
    end
    
    [temp2,Tq2] = shift(cp.supineq(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l2dc);
    Tq = Tq2(1:l2dc,:);
    
    clear temp2 Tq2
    
    C = [T(1:ldc,1:l2dc);              % Diagonal constaints
        T(ld+1:ld+k,1:l2dc)];   % Off-diagonal constaints.
    
    rel.ycons = [rel.ycons,C*temp >= 0]; 
    
    rel.A{end+1} = -(C*Tq)';
    
    clear temp Tq C dc ldc l2dc 
end

end

function rel = sdd(cp,rel) % SDD constraints.

% Declare shorthands.

n = cp.nvar; 
d = rel.relorder;
y = rel.yvar;

ld = nchoosek(n+d,d);
l2d = nchoosek(n+2*d,2*d);
tab = ncktab(n+2*d);

% General moment constraint.

mon = zeros(n,l2d);   %l2d instead of ld here because we use mon to construct the localising matrices later on.
r = zeros(n,1);

for i = 1:ld
    mon(:,i) = grlext(n,i,tab);
    r(i) = igrlext(2*mon(:,i),tab);
end

l = 0; Ay = []; Tt{(ld+1)*ld/2-ld} = [];
for i = 1:nchoosek(n+d,d)
    for j = 1:i-1
        l = l + 1;
        
        k = igrlext(mon(:,i)+mon(:,j),tab);
        
        T = zeros(3,l2d);
        T(1,r(i)) = 1; T(1,r(j)) = 1;
        T(2,k) = 2;
        T(3,r(i)) = 1; T(3,r(j)) = -1;
        
        Ay = [Ay,T*y];      % Ask Johan Lofber whether there's anything we can do to speed this up.
        Tt{l} = -T';
    end
end

rel.ycons = [rel.ycons,cone(Ay)];

rel.A{1} = Tt;  

clear T Tt Ay

% Localising matrices constraints.
    
for i = 1:numel(cp.supineq)
    
    % Declar shorthands.
    
    dc = deg(cp.supineq(i));
    ldc = nchoosek(n+floor(d-dc/2),floor(d-dc/2)); 
    l2dc = nchoosek(n+2*floor(d-dc/2),2*floor(d-dc/2)); %l2dn
    
    [temp2,Tq2] = shift(cp.supineq(i),y);   % Shift the sequence by the support polynomial.
    temp = temp2(1:l2dc);
    Tq = Tq2(1:l2dc,:);
    
    l = 1; Ay = [];
    
    for j = 1:ldc
        for k = 1:j-1

            I = igrlext(mon(:,j)+mon(:,k),tab);
            
            T = zeros(3,l2dc);
            T(1,r(j)) = 1; T(1,r(k)) = 1;
            T(2,I) = 2;
            T(3,r(j)) = 1; T(3,r(k)) = -1;
            
            Ay = [Ay,T*temp];
            Tt{l} = -Tq'*T';
            
            l = l + 1;
        end
    end
    
    rel.ycons = [rel.ycons,cone(Ay)];
    rel.A{end+1} = Tt;  
    clear T Tt Ay temp dc ldc l2dc
end

end

function rel = fwk(rel,cp)

% THIS NEEDS FIXING
% 
% % Declare shorthands.
% 
% n = cp.nvar; d = cp.relorder;
% y = cp.yvar; FW = cp.FW;
% 
% % Moment matrix constraint.
% 
% B = hankelbasis(n,d);
% 
% M = zeros(size(B{1}));
% for i = 1:nchoosek(n+2*d,2*d)
%     M = M - (-B{i})*y(i);
% end
% 
% l = nchoosek(n+d,d);
% I = nchoosek(1:l,FW);
% 
% for j = 1:numel(I(:,1))
%     cp.ycons = [cp.ycons,M(I(j,:),I(j,:)) >= 0];
% end
% 
% clear temp Ay l I
% 
% %% THIS SHOULD WORK (AND REDUCE THE SPEED ISSUES) BUT YALMIP RUNS OUT OF MEMORY.
% 
% % Ays = []; X = []; Y = [];
% % for i = 1:numel(I(:,1))
% %     Ays = [Ays;vec(M(I(i,:),I(i,:)))];
% %     X = [X;(i-1)*FW*ones(FW^2,1)+vec(repmat((1:FW)',[1,FW]))];
% %     Y = [Y;(i-1)*FW*ones(FW^2,1)+vec(repmat(1:FW,[FW,1]))];
% % end
% % 
% % cp.ycons = [cp.ycons,sparse(X,Y,Ays,numel(I(:,1))*5,numel(I(:,1))*5) >= 0];
% 
% clear l I Ay X Y
% 
% % Localising matrices constraints.
% 
% for i = 1:numel(cp.supineq)
%     
%     % Construct localising matrix.
%     
%     dc = deg(cp.supineq(i));
%     B = hankelbasis(n,floor(d-dc/2));
%     [temp,T] = shift(cp.supineq(i),y);   % Shift the sequence by the support polynomial.
%     
%     % Back to constructing the localising matrix.
%     
%     Ay = zeros(size(B{1}));    
%     
%     for j = 1:nchoosek(n+2*floor(d-dc/2),2*floor(d-dc/2))
%         Ay = Ay - (-B{j})*temp(j);
%     end
% 
%     l = nchoosek(n+floor(d-dc/2),floor(d-dc/2));
%     I = nchoosek(1:l,FW);
% 
%     if l < FW
%         cp.ycons = [cp.ycons,Ay >= 0];
%     else
%         for j = 1:numel(I(:,1))
%             cp.ycons = [cp.ycons,Ay(I(j,:),I(j,:)) >= 0];
%         end
%     end
%     clear temp Ay l I dc
% end


end

function rel = psd(cp,rel) % PSD constraints.

% Declare shorthands.

n = cp.nvar; 
d = rel.relorder;
y = rel.yvar;

ld = nchoosek(n+d,d);
l2d = nchoosek(n+2*d,2*d);
tab = ncktab(n+2*d);

% Moment matrix constraint.

B = hankelbasis(n,d);

% We need this when we recover the dual residues later.

for i = 1:l2d;
    temp{i} = - B{i};
end
rel.A{1} = temp;    
clear temp

Ay = zeros(ld);
for i = 1:l2d
    Ay = Ay - (-B{i})*y(i);
end

rel.ycons = [rel.ycons,Ay >= 0]; % C - A^T(y) >= 0

clear B; clear Ay;

% Localising matrices constraints.

for i = 1:numel(cp.supineq)
    
    % Construct localising matrix.
    
    dc = deg(cp.supineq(i));
    ldc = nchoosek(n+floor(d-dc/2),floor(d-dc/2));
    l2dc = nchoosek(n+2*floor(d-dc/2),2*floor(d-dc/2));
    
    B = hankelbasis(n,floor(d-dc/2));
    [temp,T] = shift(cp.supineq(i),y);   % Shift the sequence by the support polynomial.
    
    % We need the follwing when we recover the dual residues later.
    
    for j = 1:l2d
        temp2{j} = zeros(ldc);
        for k = 1:l2dc
            Tt = T';
            temp2{j} = temp2{j} - Tt(j,k)*B{k};
        end
    end
    
    rel.A{end+1} = temp2;    
    clear temp2

    % Back to constructing the localising matrix.
    
    Ay = zeros(ldc);    
    
    for j = 1:l2dc
        Ay = Ay - (-B{j})*temp(j);
    end
    
    rel.ycons = [rel.ycons,Ay >= 0]; % Add localising matrix constraint.
    
    clear temp Ay dc ldc l2dc
end

end

