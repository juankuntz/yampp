function plotsols(varargin)

% Creates scatter plots of the primal values of relaxations of cp that
% were succesfully solved.

% Juan Kuntz, 14/03/2015, last edited 15/05/2015

p =[];
for i = 1:nargin
    if isa(varargin{i},'scp')
        cp = varargin{i};
    elseif isa(varargin{i},'scp')
        p = varargin{i};
    end
end

if isempty(p)
    p = cp.obj{2};
end

% Remove repeated objectives

temp = p; p = [];
for i = 1:numel(temp)
    flg = 1;
    for j = 1:numel(p)
        if temp(i) == p(j)
            flg = 0;
        end
    end
    
    if flg
        p = [p;temp(i)];
    end
end

clear temp;

types = {'D','DD','SDD','FWK','PSD','NN'};

for i = 1:numel(p)
    
    % Construct arrays corersponding to finite/infinite and min/max
    % solutions.
    
    present = []; fin = []; mfin = []; inf = []; minf = [];
    
    MAX = -Inf; MIN = Inf;
    for j = 1:numel(types)
        temp = getsols(cp,{types{j},p(i)});
        [tfin,tmfin,tinf,tminf] = sortsols(temp);
        
        if ~isempty(tfin) || ~isempty(tmfin) || ~isempty(tinf) || ~isempty(tminf)
            fin{end+1} = tfin; 
            mfin{end+1} = tmfin; 
            inf{end+1} = tinf; 
            minf{end+1} = tminf;  
            present{end+1} = types{j};
            if ~isempty(tfin)
                MAX = max([MAX,tfin(2,:)]);
            end
            if ~isempty(tmfin)
                MIN = min([MIN,tmfin(2,:)]);
            end
        end
        tfin = []; tmfin = []; tinf = []; tminf = [];
    end
    
    
    figure(i);
    xlabel('rord'); ylabel('pval'); title(computestring(p(i)));
    
    c =  [1,0.8,0;1,0,0;0,0.7,0;0,0,1;0,1,1];
    %c = rand(numel(present),3);
    
    hold on
    
    bf = 0.2*(MAX-MIN)/numel(present);
    for j = 1:numel(present)
        if ~isempty(fin{j})
            scatter(fin{j}(1,:),fin{j}(2,:),[],c(j,:),'o','filled'); 
        end
        if ~isempty(mfin{j})
            scatter(mfin{j}(1,:),mfin{j}(2,:),[],c(j,:),'d','filled');  
        end
        if ~isempty(inf{j})
            scatter(inf{j}(1,:),(MAX+bf*(1+numel(present)-j))*ones(1,numel(inf{j}(1,:))),[],c(j,:),'^','filled');  
        end
        if ~isempty(minf{j})
            scatter(minf{j}(1,:),(MIN-bf*(1+numel(present)-j))*ones(1,numel(minf{j}(1,:))),[],c(j,:),'v','filled');
        end
    end
    
    % Build legend, this is a hack.

    for j = 1:numel(present)
        h(j) = plot(NaN,NaN,'color',c(j,:),'linewidth',3);
    end
    
    h = [h,plot(NaN,NaN,'marker','o','color','k','linestyle','none','MarkerFaceColor','k'),plot(NaN,NaN,'marker','d','color','k','linestyle','none','MarkerFaceColor','k'),...
        plot(NaN,NaN,'marker','^','color','k','linestyle','none','MarkerFaceColor','k'),plot(NaN,NaN,'marker','v','color','k','linestyle','none','MarkerFaceColor','k')];
    hold off
    
    present{end+1} = 'sup'; present{end+1} = 'inf'; present{end+1} = 'INF'; present{end+1} = '-INF'; 
    legend(h,present{:});
end
end

function [fin,mfin,inf,minf] = sortsols(sols)

fin = []; mfin = []; inf = []; minf = [];
for i = 1:numel(sols)
    if sols{i}.info.problem == 0 
        if strcmp(sols{i}.objs,'inf') % Feasible.
            mfin(:,end+1) = [sols{i}.rord,sols{i}.pval];
        else
            fin(:,end+1) = [sols{i}.rord,sols{i}.pval];
        end
    elseif sols{i}.info.problem == 2 
        if strcmp(sols{i}.objs,'inf') % Unbounded. 
            minf(:,end+1) = [sols{i}.rord,sols{i}.pval];
        else
            inf(:,end+1) = [sols{i}.rord,sols{i}.pval];
        end
    elseif sols{i}.info.problem == 1 
        if strcmp(sols{i}.objs,'inf') % Infeasible.
            inf(:,end+1) = [sols{i}.rord,sols{i}.pval];
        else
            minf(:,end+1) = [sols{i}.rord,sols{i}.pval];
        end
    end
end
end
