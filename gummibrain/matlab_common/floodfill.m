function new_img = floodfill(img,seedpoint,newcolor)

sz = size(img);

if(ndims(img) == 3)
    imgR = floodfill(img(:,:,1),seedpoint,-1);
    imgG = floodfill(img(:,:,2),seedpoint,-1);
    imgB = floodfill(img(:,:,3),seedpoint,-1);
    
    img_mask = zeros(sz);
    img_mask(:,:,1) = imgR < 0;
    img_mask(:,:,2) = imgG < 0;
    img_mask(:,:,3) = imgB < 0;

    img_mask = all(img_mask,3);
    
    new_img = img;
    new_img(:,:,1) = new_img(:,:,1).*(~img_mask) + (img_mask)*newcolor(1);
    new_img(:,:,2) = new_img(:,:,2).*(~img_mask) + (img_mask)*newcolor(2);
    new_img(:,:,3) = new_img(:,:,3).*(~img_mask) + (img_mask)*newcolor(3);
    
    return;
end


new_img = img;

seedcolor = img(seedpoint(1),seedpoint(2));

expansion = [-1 0; 0 1; 1 0; 0 -1];
expansion_dir = [1 2 3 4];
exclude_dir = [3 4 1 2];
num_expansion = size(expansion,1);

front = seedpoint;
front_dir = 0;

while(numel(front) > 0)
    newfront = [];
    newfront_dir = [];
    
    for f = 1:size(front,1)
        
        new_img(front(f,1),front(f,2)) = -1;
        
        new_f = expansion + repmat(front(f,:),[num_expansion 1]);
        new_d = expansion_dir;
        
        %exclude points that go back in same direction we just came
        if(front_dir(f) > 0)
            ex_d = new_d == exclude_dir(front_dir(f));
            new_f = new_f(~ex_d,:);
            new_d = new_d(~ex_d);
        end

        ex_f = new_f(:,1)<1 | new_f(:,2)<1 | new_f(:,1)>sz(1) | new_f(:,2)>sz(2);

        new_f = new_f(~ex_f,:);
        new_d = new_d(~ex_f);

        new_f_idx = sub2ind(sz,new_f(:,1),new_f(:,2));

        ex_c = new_img(new_f_idx) == -1 | new_img(new_f_idx) ~= seedcolor;

        new_f = new_f(~ex_c,:);
        new_d = new_d(~ex_c);
        
        newfront = [newfront; new_f];
        newfront_dir = [newfront_dir new_d];
    end
    
    front = newfront;
    front_dir = newfront_dir;
end
