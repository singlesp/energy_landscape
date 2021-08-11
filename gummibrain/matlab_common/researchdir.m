% location of research data on each platform
function d = researchdir

switch lower(computerid)
    case {'workpc32','workpc64'}
        %d = 'c:/users/keith/research/';
        d = 'e:/research';
    case {'workpc64new'}
        d = 'd:/data';        
    case {'macbook32','macbook64'};
        d = '~/research'; 
    case {'server1','server2'}
        d = '~';
    case {'cmrr'}
        d = '~';
    otherwise
        d = '.';
end
