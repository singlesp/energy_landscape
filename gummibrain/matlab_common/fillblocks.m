function h = fillblocks(blocktimes,blockcolors,varargin)

default_params = {'linestyle','none','facealpha',.5,'hittest','off'};

p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('timescalar',1);
p.addParamValue('ylim',[]);
p.parse(varargin{:});
r = p.Results;

timescalar = r.timescalar;
yl = r.ylim;
params = struct2args(p.Unmatched);

if(any(size(blocktimes) == 1))
    blocktimes = [ones(numel(blocktimes),1) blocktimes(:)];
end

if(any(blocktimes(:,1) ~= floor(blocktimes(:,1))))
    [a,~,u] = unique(blocktimes(:,1));
    u(blocktimes(:,1) == 0) = 0;
    count=0;
    for i = 1:numel(a)
        if(a(i) == 0)
            continue;
        end
        count = count + 1;
        u(blocktimes(:,1)==a(i)) = count;
    end
    blocktimes(:,1) = u;
end

if(nargin < 2 || isempty(blockcolors))
    blockcolors = get(gca,'colororder');
end

if(isempty(yl))
    yl = get(gca,'ylim');
end

h = zeros(size(blocktimes,1)-1,1);
count = 0;
for b = 1:size(blocktimes,1)-1
    if(blocktimes(b,1) == 0)
        continue;
    end
    count = count + 1;
    h(count) = fill(blocktimes([b b+1 b+1 b],2)*timescalar,yl([1 1 2 2]),blockcolors(blocktimes(b,1),:),default_params{:},params{:});
end
h = h(1:count);
