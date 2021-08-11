function h = barpatch(varargin)
M = varargin;
elchar = cellfun(@ischar,M);
M(elchar) = regexprep(M(elchar),'^color$','facecolor','ignorecase');
M(elchar) = regexprep(M(elchar),'^alpha$','facealpha','ignorecase');
%M = strrep(M,'color','facecolor');

mstr = cellfun(@ischar,M);
alphaprop = false(1,numel(M));
alphaprop(mstr) = regexpmatch(M(mstr),'alpha');
alphalist = {};
if(~isempty(alphaprop))
    alphaname_idx = find(alphaprop);
    alphaval_idx = find(alphaprop)+1;
    
    
    alphalist = cell(1,numel(alphaname_idx)*2);
    alphalist(1:2:end) = M(alphaname_idx);
    alphalist(2:2:end) = M(alphaval_idx);
  
    other_idx = 1:numel(M);
    other_idx([alphaname_idx alphaval_idx]) = [];    
    M = M(other_idx);
end

x=M{1};
y=M{2};

if(numel(M)<3 || isempty(M{3}))
    w=1;
else
    w=M{3};
end

if(w<1)
    h=bar(M{:});
    if(~isempty(alphalist))
        hp = findobj(h,'type','patch');
        set(hp,alphalist{:});
    end

else
    x=[x; x];
    x=x(2:end);
    y=[y;y];
    y=y(1:end-1);
    y(end)=0;

    x=[x x(end:-1:1)];
    y=[y zeros(size(y))];

    x=x(:);
    y=y(:);
    %hl=plot(x,y);
    %set(hl,'HandleVisibility','off');

    h=patch(x,y,ones(size(x)),'linestyle','-');
    set(h,M{4:end},alphalist{:});
end