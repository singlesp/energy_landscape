function cc = corrfilter2(A, template)
% normalized template match for images
% from Yue Wu, Tufts Univ ECE

mA = conv2(A,ones(size(template))./numel(template),'same');
mT = mean(template(:));
c1 = conv2(A,fliplr(flipud(template-mT)),'same')./numel(template);
c2 = mA.*sum(template(:)-mT);
stdA = sqrt(conv2(A.^2,ones(size(template))./numel(template),'same')-mA.^2);
stdT = std(template(:));
cc = (c1-c2)./(stdA.*stdT);