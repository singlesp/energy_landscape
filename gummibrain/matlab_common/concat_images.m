function Fcat = concat_images(Fimg,bgcolor,padding,direction,alignment)

if(nargin < 2 || isempty(bgcolor))
    bgcolor = 0;
end
if(nargin < 3 || isempty(padding))
    padding = 0;
end
if(nargin < 4 || isempty(direction))
    direction = 'horz';
end
if(nargin < 5 || isempty(alignment))
    alignment = 5;
end

if(strncmp(direction,'horz',4))
    dir_horz = true;
else
    dir_horz = false;
end

Fimg = cellfun(@(x)(padimage(x,padding,bgcolor)),Fimg,'uniformoutput',false);
maxsz = max(cell2mat(cellfun(@size,Fimg,'uniformoutput',false)'),[],1);

for i = 1:numel(Fimg)
    imgsz = size(Fimg{i});
    if(dir_horz)
        imgsz(1) = maxsz(1);
    else
        imgsz(2) = maxsz(2);
    end
    Fimg{i} = padimageto(Fimg{i},imgsz,alignment,bgcolor);
end

if(dir_horz)
    Fcat = [Fimg{:}];
else
    Fimg = cellfun(@(x)(permute(x,[2 1 3])),Fimg,'uniformoutput',false);
    Fcat = [Fimg{:}];
    Fcat = permute(Fcat,[2 1 3]);
end
