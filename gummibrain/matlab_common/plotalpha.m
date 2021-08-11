function h = plotalpha(varargin)
% works like normal plot except it uses fill to allow alpha transparency
%
% plotalpha(x,y,'r')
% plotalpha(x,y,'alpha',.5)


args = {};

if(nargin >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2}))
    x = varargin{1};
    y = varargin{2};
    args = varargin(3:end);
elseif(nargin >= 2 && isnumeric(varargin{1}) && ~isnumeric(varargin{2}))
    %x = varargin{1};
    x = [];
    y = varargin{1};
    args = varargin(2:end);
elseif(nargin == 1 && isnumeric(varargin{1}))
    x = [];
    y = varargin{1};
    args = varargin(2:end);
else
    error('invalid inputs');
end

def_edgecolor = 'r';
def_alpha = .15;


if(~isempty(args) && ischar(args{1}) && numel(args{1}) == 1)
    def_edgecolor = args{1};
    args = args(2:end);
end
% 
% if(~isempty(args) && isnumeric(args{1}) && numel(args{1} == 3))
%     def_edgecolor = args{1};
%     args = args(2:end);
% end

charidx = cellfun(@ischar,args);
args(charidx) = regexprep(args(charidx),'^color$','edgecolor');
args(charidx) = regexprep(args(charidx),'^alpha$','edgealpha');

do_markers = false;
markeralpha = def_alpha;
markeridx = strcmpi(args(charidx),'markeralpha');
if(any(markeridx))
    ch = element(find(charidx),find(markeridx,1,'first'));
    markeralpha = args{ch+1};
    if(markeralpha > 0)
        do_markers = true;
    end
    args = args(setdiff(1:numel(args),[ch ch+1]));
end

if(size(y,1) == 1)
    y = reshape(y,[],1);
end
if(isempty(x))
    x = [1:size(y,1)].';
end

if(size(x,1) == 1)
    x = reshape(x,[],1);
end

if(numel(x) == size(y,1))
    x = repmat(x,1,size(y,2));
end

lx = [x; x(end,:)+1];
ly = [y; nan(1,size(y,2))];

h=fill(lx(:),ly(:),'w','facealpha',0,'edgecolor',def_edgecolor,'edgealpha',.15,'linewidth',1,args{:});

if(do_markers)
    circN = 10;
    circR = .01;
    cx = circR*cos(linspace(0,2*pi,circN));
    cy = circR*sin(linspace(0,2*pi,circN));
    
    mx = repmat(cx(:),1,numel(x)) + repmat(reshape(x,1,[]),circN,1);
    my = repmat(cy(:),1,numel(y)) + repmat(reshape(y,1,[]),circN,1);
    

    np = get(gca,'nextplot');
    hold on;
    h2=fill(mx,my,def_edgecolor,'linestyle','none','facealpha',markeralpha);
    set(gca,'nextplot',np);    
    h = [h; h2];
end
