function varargout = complex_ellipse(cplxdata,varargin)

p = inputParser;
p.addParamValue('color',[1 0 0]);
p.addParamValue('alpha',.95);
p.addParamValue('std',[]);
p.addParamValue('plotsamples',false);
p.addParamValue('blocksize',1);
p.addParamValue('fillstyle',{});
p.addParamValue('plotstyle',{});
p.addParamValue('draw',true);

p.parse(varargin{:});
r = p.Results;

color = r.color;
do_plot = r.plotsamples;
blocksize = r.blocksize;
plotstyle = r.plotstyle;
fillstyle = r.fillstyle;

do_draw = r.draw;

if(~isempty(r.std))
    ci_alpha = 1-(normcdf(-r.std,0,1)*2);
else
    ci_alpha = r.alpha;
end

linestyle = {};
pointstyle = {'marker','.'};

p=2; %dimensionality
n=size(cplxdata,1);
k=finv(ci_alpha,p,n-p)*p*(n-1)/(n-p);

ellipseinfo = {};

for i = 1:size(cplxdata,2)
    data1 = [real(cplxdata(:,i)) imag(cplxdata(:,i))];
    ellipseinfo{i} = confidence_ellipsoid(data1,'alpha',ci_alpha);
end

if(do_draw)
    axmax = 0;
    set(gca,'nextplot','add');
    for i = 1:numel(ellipseinfo)
        for j = 1:numel(ellipseinfo{i});
            m1 = ellipseinfo{i}(j).mean;
            el1 = ellipseinfo{i}(j).ellipsexy;
            fill(m1(1)+el1(:,1),m1(2)+el1(:,2),color,'linestyle','none','facealpha',.5,fillstyle{:});            
        end

        [~,midx] = max([ellipseinfo{i}.alpha]);
        m1 = ellipseinfo{i}(midx).mean;
        ax1 = ellipseinfo{i}(midx).axpc;
        el1 = ellipseinfo{i}(midx).ellipsexy;
        
        m = max(flatten(abs([repmat(m1,size(el1,1),1)+el1])));
        axmax = max(axmax,m);
        
        plot([0; m1(1)], [0; m1(2)],'color',color);
        plot([-ax1(:,1) ax1(:,1)]'+m1(1), [-ax1(:,2) ax1(:,2)]'+m1(2),'k');
        

        if(do_plot)
            dsz = blocksize*floor(size(cplxdata,1)/blocksize);
            blockdata = reshape(cplxdata(1:dsz,i),[],blocksize);
            blockdata = mean(blockdata,2);
            plot([zeros(size(blockdata,1),1) blockdata].','color',color, linestyle{:});
            plot(blockdata,'color',color, pointstyle{:});
        end
    end
    
    %axis tight;
    %m = max(flatten(abs(Ftmp)));
    %m = max(abs([get(gca,'xlim') get(gca,'ylim')]));
    %m = max(flatten(abs([repmat(m1,size(el1_outer,1),1)+el1_outer])));
    if(axmax > 0)
        axis([-axmax axmax -axmax axmax]);
    end
    axis square;
    %set(gca,'xtick',[0],'ytick',[0],'xticklabel',[],'yticklabel',[]);
    grid on;
    %return;
    %axis tight;    
end

if(nargout > 0)
    varargout = ellipseinfo;
end




