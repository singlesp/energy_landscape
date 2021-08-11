function ax = plotradial(datacplx,pointstyle)

markertypes = '+o*.xsd^v><ph';
colortypes = 'rgbcmykw';
linetypes = {'-', '--', ':', '-.'};

 pointstyle = {'g-.','color'};
    


ltype = '';
mtype = '';
mcolor = '';
if(~isempty(pointstyle) && ischar(pointstyle{1}))
    lfind = cellfun(@(x)(strfind(pointstyle{1},x)),linetypes,'uniformoutput',false);
    lidx = find(~cellfun(@isempty,lfind),1,'last');
    if(~isempty(lidx))
        ltype = linetypes{lidx};
        pointstyle{1}(lfind{lidx} + (0:numel(ltype)-1)) = '';
    end
    mtype_tmp = intersect(pointstyle{1},markertypes);
    mcolor_tmp = intersect(pointstyle{1},colortypes);
    if(~isempty(mtype_tmp))
        mtype = mtype_tmp(1);
    end
    if(~isempty(mcolor_tmp))
        mcolor = mcolor_tmp(1);
    end
    pointstyle{1}(pointstyle{1} == mtype | pointstyle{1} == mcolor) = '';
    if(isempty(pointstyle{1}))
        pointstyle = pointstyle(2:end);
    end
end

if(isempty(mtype))
    mtype = 'none';
end
if(isempty(ltype))
    ltype = 'none';
end

pointstr = {'markertype',mtype,'color',mcolor};
if(isempty(ltype))
    linestr = {'linestyle','none'};
else
    linestr = {'linestyle',ltype,'color',mcolor};
end

%co = [1 0 0; 0 0 1];
co = get(gca,'colororder');
hold on;
for i = 1:size(datacplx,2)
    if(isempty(mcolor))
        c = co(i,:);
    else
        c = mcolor
    
    plot(datacplx(:,i),'.','color',c,pointstr{:});

    plot([zeros(size(datacplx,1),1) datacplx(:,i)].','color',c,'markerstyle','none',pointstyle{:});
end

m = max(abs(datacplx(:)));
axis([-m m -m m]);
axis square;
set(gca,'xtick',[0],'ytick',[0]);
grid on;