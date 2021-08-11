function yvals=crosspoint(ax,xvals,varargin)

gridmode=false;
gridmode_idx=find(cellfun(@(x)(strcmpi(x,'gridmode')),varargin));

if(~isempty(gridmode_idx))
    newidx=1:numel(varargin);
    rmidx=gridmode_idx;
    if(numel(varargin) > gridmode_idx)
        if(strcmpi(varargin{gridmode_idx+1},'on'))
            gridmode=true;
        end
        rmidx=[rmidx gridmode_idx+1];
    end
    newidx(rmidx)=[];
    varargin=varargin(newidx);
end

for a = ax
    np=get(a,'nextplot');
    set(a,'nextplot','add');
    
    xl = get(a,'xlim');
    yl = get(a,'ylim');
    
    for v = 1:numel(xvals)
        
        xval=xvals(v);
        yvals=getcrosspoint(a,xval);
        for i = 1:numel(yvals)
            yval=yvals{i};
            if(gridmode && i==1)
                plot(a,[xval xval],[yl(1) yl(2)],varargin{:},'handlevisibility','off');
            end
            if(~isempty(xval) && ~isempty(yval))
                if(gridmode)
                    plot(a,[xl(1) xl(2)],[yval yval],varargin{:},'handlevisibility','off');
                else
                    plot(a,[xval xval],[yl(1) yval],varargin{:},'handlevisibility','off');
                    plot(a,[xl(1) xval],[yval yval],varargin{:},'handlevisibility','off');
                end
            end
        end
    end
    
    set(a,'nextplot',np);
end