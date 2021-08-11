function compstr = computerid

archstr = upper(computer);
compname = upper(getenv('COMPUTERNAME'));
hostname = upper(getenv('HOSTNAME'));
if(isempty(hostname))
  [~,hostname] = system('hostname');
end

compstr = ['other_' archstr];

%some known computers
if(strcmp(archstr,'PCWIN') && strcmp(compname,'KEITH-PC'))
    compstr = 'workpc32';
    
elseif(strcmp(archstr,'PCWIN64') && strcmp(compname,'KEITH-PC'))
    compstr = 'workpc64';
    
elseif(strcmp(archstr,'PCWIN64') && strcmp(compname,'BME-HE6112KJ2'))
    compstr = 'workpc64new';
    
elseif(strcmp(archstr,'MACI'))
    compstr = 'macbook32';
    
elseif(strcmp(archstr,'MACI64'))
    compstr = 'macbook64';
    
elseif(strcmp(archstr,'GLNXA64') && ~isempty(regexpi(hostname,'^han_?server1')))
    compstr = 'server1';
    
elseif(strcmp(archstr,'GLNXA64') && ~isempty(regexpi(hostname,'^han_?server2')))
    compstr = 'server2';
    
elseif(strcmp(archstr,'GLNXA64') && ~isempty(regexpi(hostname,'^.+.cmrr.umn.edu')))
    compstr = 'cmrr';
end
