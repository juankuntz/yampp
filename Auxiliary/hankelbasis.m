function B = hankelbasis(n,d)

% B{:} is the "moment" basis for the space of matrices with the same type
% of Hankel-like structure that moment matrices (of dimension n+d choose d)
% have.

% Juan Kuntz, 12/02/2015

nd = nchoosek(n+d,d); n2d = nchoosek(n+2*d,2*d);

B{n2d} = [];

inds = summi(n,d);

B{n2d} = [];

for i = 1:n2d
    B{i} = sparse(inds{i}(1,:),inds{i}(2,:),ones(1,numel(inds{i}(1,:))),nd,nd,numel(inds{i}(1,:)));
end
end

