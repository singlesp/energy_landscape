function imgblocks = BreakImageIntoBlocks(img, blocksize)

[imgw imgh] = size(img);

imgblocks = {};

xstart = 1:blocksize:imgw;
ystart = 1:blocksize:imgh;

for i = 1:numel(xstart)-1
    for j = 1:numel(ystart)-1
        imgblocks{numel(imgblocks)+1} = img(xstart(i):xstart(i)+blocksize-1,ystart(i):ystart(i)+blocksize-1);
    end
end