% new_img = CropBGColor(img, bgcolor)
% bgcolor 
function [new_img croprect] = CropBGColor(imgs, bgcolor)

if(~iscell(imgs))
    imgs = {imgs};
end
    
new_img = {};
for i = 1:numel(imgs)
    img = imgs{i};
    
    nd = ndims(img);
    if(~(nd == 2 || nd == 3))
        error('must be 2 or 3 dimensional');
    end

    if(nargin == 1)
        if(nd == 2)
            bgcolor = 255;
        elseif(nd == 3)
            bgcolor = [255 255 255];
        end
    end

    bgcolor = squeeze(bgcolor);
    if(numel(bgcolor) < nd)
        bgcolor = repmat(bgcolor,1,nd);
    end

    if(nd == 2)
        bg = (img == bgcolor);
        [fgx fgy] = ind2sub(size(bg),find(~bg));
        new_img{i} = img(min(fgx):max(fgx), min(fgy):max(fgy));
    elseif(nd == 3)
        bg = (img(:,:,1) == bgcolor(1)) & (img(:,:,2) == bgcolor(2)) & (img(:,:,3) == bgcolor(3));
        [fgx fgy] = ind2sub(size(bg),find(~bg));
        new_img{i} = img(min(fgx):max(fgx), min(fgy):max(fgy),:);
    end
end
if(numel(new_img) == 1)
    new_img = new_img{1};
end

if(nargout > 1 && numel(imgs) == 1)
    croprect = [min(fgx) min(fgy) max(fgx) max(fgy)];
end