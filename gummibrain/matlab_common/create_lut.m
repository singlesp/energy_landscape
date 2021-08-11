lut = [255 0 0; 0 0 0; 0 255 0];
N = size(lut,1);
lutB = interp1(1:N,lut,linspace(1,N,256));
lutB = [[0:size(lutB,1)-1]' lutB];

fprintf('* s=byte	index	red	green	blue   alpha\n');
fprintf('S\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n',lutB')
