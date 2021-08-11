function output_txt = datatip_showindex(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

target = get(event_obj,'Target');
type = lower(get(target,'Type'));
dataidx = nan;
dataval = nan;

switch(type)
    case {'patch'}
         verts = get(target,'Vertices');
         vals = get(target,'FaceVertexCData');
         [dataidx err] = FindClosestPoints(pos,verts);
         if(numel(vals) == size(verts,1))
             dataval = vals(dataidx);
         end
%         disp(['Index: ' num2str(dataidx)]);
         
%          %apparently can't do this without triggering extra datatips         
%          ax = get(target,'Parent');
%          fig = get(ax,'Parent');
%          nextplot = get(ax,'NextPlot');
%          set(ax,'NextPlot','add');
%          plot3(verts(dataidx,1),verts(dataidx,2),verts(dataidx,3),'ko',...
%             'MarkerFaceColor','w','MarkerSize',10);
%          set(ax,'NextPlot',nextplot);
    case {'surface'}
        c = get(target,'CData');
         x = get(target,'XData');
         y = get(target,'YData');
         z = get(target,'ZData');
         verts = [x(:) y(:) z(:)];
         [dataidx err] = FindClosestPoints(pos,verts);
         %dataidx = c(dataidx);
         dataval = c(dataidx);
         %disp(['Index: ' num2str(c(dataidx))]);
    case 'line'
        dataidx = get(event_obj,'DataIndex');
end



output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

if(~isnan(dataidx))
    output_txt{end+1} = ['Index: ', num2str(dataidx)];
end

if(~isnan(dataval))
    output_txt{end+1} = ['Value: ', num2str(dataval)];
end
