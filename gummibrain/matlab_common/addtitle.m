%new_img = addtitle(img,title,fontsize,bgcolor)
function new_img = addtitle(img,titlestr,fontsize,bgcolor)

if(nargin < 3)
    fontsize = 20;
end

if(nargin < 4)
    if(ndims(img)) == 3
        bgcolor = [255 255 255];
    else
        bgcolor = 255;
    end
end

new_img = padimage(img, [fontsize + 10 0 0 0], bgcolor);
new_img = addtext2img(new_img, {{size(new_img,2)/2, size(new_img,1), ...
        titlestr, ...
       'VerticalAlignment','top','HorizontalAlignment','center','Units','Pixels','FontUnits','pixels','FontSize',fontsize}});