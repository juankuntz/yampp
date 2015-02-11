limpio

n = 6; d = 4;

tic 
for k = 1:100
for i = 1:50
    for j = 1:i
        nchoosek(i,j);
    end
end
end
toc

tic 
for k = 1:100
for i = 1:50
    for j = 1:i
        binomial(i,j);
    end
end
end

toc
% A = [];
% tic 
% for e = 0:d
%     A = [A,genmulti(n,e)];
% end
% toc
% 
% B = [];
% tic 
% for e = 0:d
%     B = [B,grlexb(n,e)];
% end
% toc

% A = []; B = [];
% 
% tic
% for i = 1:nchoosek(n+d,d)
%     A = [A,grlex(n,i)];
% end
% toc
% 
% tic
% for i = 1:nchoosek(n+d,d)
%     B = [B,grlex(n,i)];
% end
% toc

% pol('x',n); pol('y',n);
% mpol('X',n); mpol('Y',n);
% 
% 
% cx= rand(1,n); cy = rand(1,n);
% 
% tic 
% a = (cx*X+cy*Y)^d;
% toc
% 
% tic 
% b = (cx*x+cy*y)^d
% toc
% 
