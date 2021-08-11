function hfig = resizefig(hfig,vsize,hsize)
% resizefig(gcf,2) %scale figure
% resizefig(gcf,2,1)
% resizefig(gcf,'screen') %scale relative to screen size
% resizefig(gcf,'screen*.5')
% resizefig('all','screen*.5') %all figures


if(~any(ishandle(hfig)) && ischar(hfig))
    hfig = findobj(0,'type','figure');
end

for i = 1:numel(hfig)
    if(~ishandle(hfig(i)))
        continue;
        
    end
    if(nargin == 2)
        if(isnumeric(vsize) && numel(vsize) == 2)
            hsize = vsize(2);
            vsize = vsize(1);
        else
            hsize = vsize;
        end
    end
    
    u = get(0,'units');
    set(0,'units','pixels');
    screensz = get(0,'ScreenSize');
    set(0,'units',u);
    figsz = getpixelposition(hfig(i));
    
    if(ischar(vsize))
        vsize = eval(regexprep(vsize,'screen',sprintf('%d',screensz(4))));
    end
    
    if(isnumeric(vsize))
        vsize = figsz(4)*vsize;
    end
    
    if(ischar(hsize))
        hsize = eval(regexprep(hsize,'screen',sprintf('%d',screensz(3))));
    end
    if(isnumeric(hsize))
        hsize = figsz(3)*hsize;
    end
    
    newsz = [figsz(1:2) hsize vsize];
    setpixelposition(hfig(i),newsz);
end