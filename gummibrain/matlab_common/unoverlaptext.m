function unoverlaptext(varargin)
va = varargin;
offsetamt = 0;
idx = StringIndex(va,'offset');
if(numel(va) > 0 && ~isempty(idx))
    offsetamt = va{idx+1};
    va = va(setdiff(1:numel(va),[idx idx+1]));
end

if(numel(va) > 0 && all(strcmpi(get(va{1},'type'),'axes')))
    ax = va{1};
elseif(numel(va) > 0 && all(strcmpi(get(va{1},'type'),'figure')))
    ax = gaa(va{1});
end

for a = 1:numel(ax)
    ht = findobj(ax(a),'type','text');
    np = get(ax(a),'nextplot');
    u = get(ax(a),'units');
    set(ax(a),'nextplot','add');
    %set(ax(a),'units','pixels');
    
    if(numel(ht) == 0)
        continue;
    elseif(numel(ht) == 1)
        pos = get(ht,'position');
    else
        pos = cell2mat(get(ht,'position'));
    end
    [~,sortidx] = sort(pos(:,2),'ascend');
    ht = ht(sortidx);
    
    r = [];
    for i = 1:numel(ht)
        ut = get(ht(i),'units');
        %set(ht(i),'units','pixels');
        pos = get(ht(i),'position');
        ext = get(ht(i),'extent');
        ext(3:4) = ext(3:4)+ext(1:2);
        ext_orig = ext;
        if(size(r,1) > 0)
            o1 = (ext(1) >= r(:,1) & ext(1) <= r(:,3));
            o2 = (ext(3) >= r(:,1) & ext(3) <= r(:,3));
            o3 = (ext(2) >= r(:,2) & ext(2) <= r(:,4));
            o4 = (ext(4) >= r(:,2) & ext(4) <= r(:,4));

            o5 = (ext(1) <= r(:,1) & ext(3) >= r(:,3));
            o6 = (ext(2) <= r(:,2) & ext(4) >= r(:,4));
            
            o1 = o1 | o5;
            o2 = o2 | o5;
            o3 = o3 | o6;
            o4 = o4 | o6;
            o = (o1 & o3) | (o1 & o4) | (o2 & o3) | (o2 & o4);
            
            %o = ext(1) <= r(:,3) & ext(3) >= r(:,1) & ext(2) <= r(:,4) & ext(4) >= r(:,2);
            
            if(any(o))
                dy = max(r(o,4))-ext(2);
                ext([2 4]) = ext([2 4]) + dy;
            end
        end

        newpos = pos;
        newpos(1) = newpos(1) + offsetamt;
        if(~all(ext==ext_orig))
            %newpos(1) = newpos(1) + offsetamt;
            newpos(2) = newpos(2) + (ext(2)-ext_orig(2));
        end
        
        set(ht(i),'position',newpos);
        
        plot(ax(a),pos(1) + [0 .75*(newpos(1)-pos(1))],[pos(2) newpos(2)],'color',get(ht(i),'color'));
        r = [r; ext];
        set(ht(i),'units',ut);
    end
    set(ax(a),'nextplot',np,'units',u);
end