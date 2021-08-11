% location of matlab code on each platform
function d = codedir

switch lower(computerid)
    case {'workpc32','workpc64'}
        d = 'c:/users/keith/research/code';
    case {'macbook32','macbook64'};
        d = '~/research/code'; 
    case {'server1','server2'}
        d = '~/code';
    otherwise
        d = './';
end
