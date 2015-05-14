function varargout = shift(p,y)

% Shift sequence varargin{2} (of class seq or sdpvar) by polynomial varargin{1}.

% Currently it automatically matches the first variable of p with that of y
% and so on. 

% Warning: Because of how the sdpvar class currently works, if the output
% of this function is not stored in a variable whose name is a string of
% length greater than one, then sdpvar will start remaining variables which
% may confuse things. For example, if we pass out a single output y(2) and
% this is stored in z = shift(p,y), then y(2) in y will be renamed as z...
%
% Edit: Has the option of returning a second argument T that is the
% transformation matrix from y to z = shift(p,y), that is z = Ty. This is
% useful when computing the problem data.
%
% This needs cleaning
%
% Juan Kuntz, 13/02/2015, last edited 14/05/2015.

% Check for valid arguments.

if (~isa(y,'sdpvar') && ~isa(y,'seq')) || ~isa(p,'pol')
    disp('Error: a pol object and either a seq object r passed when calling /yampp/@pol/shift');
    return
elseif isa(y,'seq') && y.dim ~= p.nvar
    disp('Error: The dimension of the underlying space of the sequence must be the same as that of the polynomial');
    return
end

if p == pol(1)
    varargout{1} = y;
    varargout{2} = eye(numel(y));
    return
end

if isa(y,'sdpvar') % Compute the order of sequence.
    d = sum(grlex(p.nvar,numel(y)));
else
    d = y.ord;
end

T = zeros(nchoosek(p.nvar+d-p.deg,d-p.deg),nchoosek(p.nvar+d,d)); % Initialise

% Compute the exponent vectors corresponding to every monomial in the
% support of p.

tab = ncktab(p.nvar+p.deg+d);  % Compute appropiate table for fast grlex igrlex use.

temp = zeros(p.nvar,numel(p.coef(1,:)));
for i = 1:numel(p.coef(1,:))
    temp(:,i) = grlext(p.nvar,p.coef(2,i),tab);
end

% If y is sdpvar object.

if isa(y,'sdpvar')

    % out_a = sum_b (p_b)*(y_(a+b)) where p_b denote the coefficient of p
    % corresponding to monomial b.
    out = sdpvar(nchoosek(p.nvar+(d-p.deg),(d-p.deg)),1);
    
    for i=1:nchoosek(p.nvar+(d-p.deg),(d-p.deg));
        tempi = grlext(p.nvar,i,tab);   
        test = 0;
        for j = 1:numel(p.coef(1,:))
            temp2 = igrlext(tempi+temp(:,j),tab);   
            if test == 0
                test = 1;
                out(i) =  p.coef(1,j)*y(temp2);
                T(i,temp2) = T(i,temp2) + p.coef(1,j);
            else
                out(i) =  p.coef(1,j)*y(temp2) + out(i);
                T(i,temp2) = T(i,temp2) + p.coef(1,j);
            end
        end
        clear tempi
    end
    varargout{1} = out;
    varargout{2} = T;
    return
end

% If y is a seq object.

out = seq(p.nvar,d-p.deg); % Initialise seq of appropiate dimension.

ARRAY = y.coef; % The variables ARRAY and ARRAY2 are a hack needed because subsref and subsasgn for seq objects is not well developed (for example, it cannot handle y.coef(1,1) if y is seq).
ARRAY2 = zeros(2,nchoosek(p.nvar+(d-p.deg),(d-p.deg)));

% out_a = sum_b (p_b)*(y_(a+b)) where p_b denote the coefficient of p
% corresponding to monomial b.

for i = 1:nchoosek(p.nvar+(d-p.deg),(d-p.deg));
    tempi = grlext(p.nvar,i,tab);
    test = 0;
    
    for j = 1:numel(p.coef(1,:))
        temp2 = igrlext(tempi+temp(:,j),tab);
        I = bfind(ARRAY(2,:),temp2);
        if ~isempty(I)
            if test == 0
                test = 1;
                ARRAY2(:,i) = [p.coef(1,j)*ARRAY(1,I);i];
                T(i,temp2) = T(i,temp2) + p.coef(1,j);
            else
                ARRAY2(:,i) = ARRAY2(1,end) + p.coef(1,j)*ARRAY(1,I);
                T(i,temp2) = T(i,temp2) + p.coef(1,j);
            end
        end
    end
    
out.coef = ARRAY2;

varargout{1} = out;
varargout{2} = T;
end