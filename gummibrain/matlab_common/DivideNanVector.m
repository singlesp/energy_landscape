function vlim = DivideNanVector(v)

vbreaks = find(isnan(v));
if(size(vbreaks,1) > size(vbreaks,2))
    vbreaks = vbreaks';
end

if(isempty(vbreaks))
    vlim = [1 numel(v)];
else
    vlim = [1 vbreaks(1)-1; vbreaks(1:end-1)'+1 vbreaks(2:end)'-1];
end

if(~isempty(vbreaks) && vbreaks(end) < numel(v))
    vlim = [vlim; vbreaks(end)+1 numel(v)];
end