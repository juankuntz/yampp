function mkscp(cp)

% Construct the appropiate moment problem.

% Juan Kuntz, 13/02/2015

n = cp.nvar; d = cp.relorder;

% Initialise yalmip variables.

cp.yvar = [];
cp.yvar = sdpvar(nchoosek(n+2*d,2*d),1); 
y = cp.yvar;

% Clear any previous yalmip constraints.

cp.ycons = [];

% PSD constraints

% Moment matrix constraint.

B = hankelbasis(n,d);

Ay = zeros(size(B{1}));
for i = 1:nchoosek(n+2*d,2*d)
    Ay = Ay + B{i}*y(i);
end

cp.ycons = [cp.ycons,Ay >= 0];

clear B; clear Ay;

% Localising matrices constraints.

for i = 1:numel(cp.supcon)
    if 2*d < cp.supcon(i).deg
        disp('Error: degree of one of the polynomials defining the support is too big'); % fix this.
        return
    end
    
    % Construct localising matrix.
    
    B = hankelbasis(n,floor(d-cp.supcon(i).deg/2));
    
    Ay = zeros(size(B{1}));
    temp = shift(cp.supcon(i),y);
    
    for j = 1:nchoosek(n+2*floor(d-cp.supcon(i).deg/2),2*floor(d-cp.supcon(i).deg/2))
        Ay = Ay + B{j}*temp(j);
    end
    
    cp.ycons = [cp.ycons,Ay >= 0]; % Add localising matrix constraint.
    
    clear temp Ay
end

for i = 1:numel(cp.seqeqcon)
    cp.ycons = [cp.ycons,cp.seqeqcon(i)*y == 0];
end
end

