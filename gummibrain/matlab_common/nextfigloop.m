function nextfigloop(pause_sec)
if ~exist('pause_sec','var') || isempty(pause_sec)
    pause_sec=[];
end

while(true)
    nextfig;
    if(isempty(pause_sec))
        pause;
    else
        pause(pause_sec);
    end
end