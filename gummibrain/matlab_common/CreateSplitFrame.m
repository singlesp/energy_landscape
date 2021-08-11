function varargout = CreateSplitFrame(hParent,split_direction,numframes,framesizes,gapsize)

if(nargin < 3 || isempty(numframes))
    numframes = 2;
end

if(nargin < 4)
    framesizes = ones(1,numframes)/numframes;
end

framesizes = reshape(framesizes,1,[]);

if(numel(framesizes) < numframes)
    framesize_leftover = 1-sum(framesizes);
    numleftover = numframes - numel(framesizes);
    framesizes = [framesizes framesize_leftover*ones(1,numleftover)/numleftover];
end

if(nargin < 5)
    gapsize = .025;
end

frames = {};
splitframes = {};

if(strcmpi(split_direction,'horz'))
    d = [1 3];
elseif(strcmpi(split_direction,'vert'))
    d = [2 4];
end

curpos = 0;
for i = 1:numframes
    pos = [0 0 1 1];
    
    framebound = framesizes(i)-(i<numframes)*gapsize/2;

    pos(d) = [curpos framebound];
    pos(d(2)) = min(pos(d(2)),1-pos(d(1)));
    frames{i} = uipanel(hParent,'units','normalized','position',pos);
    curpos = sum(pos(d));
    if(i < numframes)
        pos = [0 0 1 1];
        pos(d) = [curpos gapsize];
        splitframes{i} = uipanel(hParent,'units','normalized','position',pos,'bordertype','none');
    end
    curpos = sum(pos(d));
end

for i = 1:numframes-1
    set(splitframes{i},'buttondownfcn',{@frameResize_click, [split_direction frames(i:i+1) splitframes{i}]},...
        'userdata',struct('pixelpos',getpixelposition(splitframes{i})));
end

if(strcmpi(get(hParent,'type'),'uipanel'))
    set(hParent,'bordertype','none');
end

hFig = GetParentFigure(hParent);
if(~isappdata(hFig,'splitframe'))
    iptaddcallback(hFig,'WindowButtonUpFcn',{@frameResize_click,[]});
    iptaddcallback(hFig,'WindowButtonMotionFcn',@frameResize_drag);
end
A = getappdata(hFig,'splitframe');
A.current_frameinfo = [];
setappdata(hFig,'splitframe',A);

if(nargout == 1)
    varargout = {frames};
elseif(nargout == 2)
    varargout = {frames, splitframes};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frameResize_click(hObj,event,frameinfo)
A = getappdata(gcbf,'splitframe');
A.current_frameinfo = frameinfo;
setappdata(gcbf,'splitframe',A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frameResize_drag(hObj,event)
A = getappdata(gcbf,'splitframe');
frameinfo = A.current_frameinfo;
if(isempty(frameinfo) || numel(frameinfo) ~= 4)
    return;
end


split_direction = frameinfo{1};
frame1 = frameinfo{2};
frame2 = frameinfo{3};
splitframe = frameinfo{4};
hParent = get(frame1,'parent');

%framedata = get(splitframe,'userdata');


split_pos = get(splitframe,'position');
frame1_pos = get(frame1,'position');
frame2_pos = get(frame2,'position');


if(strcmpi(split_direction,'horz'))
    d = [1 3];
elseif(strcmpi(split_direction,'vert'))
    d = [2 4];
end

framebounds = [frame1_pos(d(1)) sum(frame2_pos(d))];


curpos = get(hObj,'currentpoint');

if(strcmpi(get(hParent,'type'),'figure'))
    form_pos = get(hObj,'position');
    curpos = curpos ./ form_pos(3:4);
else
    container_pos = getpixelposition(hParent,true);
    curpos = curpos - container_pos(1:2);
    curpos = curpos ./ container_pos(3:4);
end


curpos(d(1)) = min(curpos(d(1)), framebounds(2)-split_pos(d(2)));
curpos(d(1)) = max(curpos(d(1)), framebounds(1)+split_pos(d(2)));

split_pos(d(1)) = curpos(d(1))-split_pos(d(2))/2;
frame2_pos(d(1)) = sum(split_pos(d));
frame2_pos(d(2)) = framebounds(2)-frame2_pos(d(1));
frame1_pos(d(2)) = split_pos(d(1))-framebounds(1);

set(splitframe,'position',split_pos);
set(frame1,'position',frame1_pos);
set(frame2,'position',frame2_pos);
