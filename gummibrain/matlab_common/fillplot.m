function h=fillplot(x,y,varargin)
args=varargin;

do_step=false;
while ~isempty(args)
    if(strcmpi(args{1},'step'))
        do_step=true;
        args=args(2:end);
    else
        break;
    end
    
end
if(size(y,1)<size(y,2))
    y=y.';
end

if(do_step)
    
    i=[1:numel(x) numel(x); 1:numel(x) 1];
    j=[1:size(y,1)-1; 2:size(y,1)];
    
    h=fill(x(i(:)),[zeros(2,size(y,2)); y(j(:),:); zeros(2,size(y,2))],args{:});
else
    h=fill(x([1:end end 1]),[y(1:end,:); zeros(2,size(y,2))],args{:});
end

