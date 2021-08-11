function varargout = unixvm(cmdstr, do_echo, dirmaps_new)
% call a linux command, locally if operating in linux, or via ssh with
%  a virtual machine.  Note: This can also be used for remote machines, but
%  the file sharing will be slightly more complicated.
%
% Before you can use this function, you need to do the following:
%
% On the remote/virtual linux machine:
%  1. Mount your local harddrives so that the remote machine and your local
%   PC are accessing the same data.  
%  2. Set the filename replacements
%  3. Run puttygen, generate ssh keypairs, save the putty key file (eg:
%  id_rsa.ppk), the openssh private key (
%  3. Run pageant (putty ssh agent.. best to have this run on Windows
%   startup)
% .....

%get the unix case out of the way...
if(~ispc)
    %status = exit code.  result = output text
    %TODO: need to find out how to get exit code from sshconn
    [status,result] = system(cmdstr);
    exitcode = status;
    result = reshape(regexp(result,'[\n]','split'),[],1);
    if(nargout == 1)
        varargout = {result};
    elseif(nargout == 2)
        varargout = {result,result};
    elseif(nargout == 3)
        varargout = {result,result,exitcode};
    end
    return;
end

%%
%remote IP, remote user, local pubkey, remote port
connparams = {'192.168.56.1','brain','C:/Users/kjamison/.ssh/id_rsa','',2222};
%%
global vmconn;
persistent dirmaps;
%%
if(exist('dirmaps_new','var'))
    dirmaps = dirmaps_new;
end

if(isempty(cmdstr))
    %called unixvm([],[],dirmap) so initialize dirmap
    varargout = {};
    return;
end

if(~exist('dirmaps','var') || isempty(dirmaps))
    %see unixvmpath_defaults.m to set the windows->linux (virtualbox)
    %mapping
    dirmaps = unixvmpath_defaults;
end


if(~exist('no_echo','var') || isempty(do_echo))
    do_echo = true;
end


%%

[args, dequoteargs] = parseargs(cmdstr);
for i = 1:numel(args)
    do_fixslash = false;
    if(any(args{i}=='\'))
        %if this is an existing file or folder
        if(any(exist(dequoteargs{i}) == [2 7]))
            do_fixslash = true;
        else
            [tmpd,~,~] = fileparts(dequoteargs{i});
            if(any(exist(tmpd) == [2 7]))
                do_fixslash = true;
            end
        end
    end
    if(do_fixslash)
        args{i} = strrep(args{i},'\','/');
    end
end
cmdstr = sprintf('%s ',args{:});

cmdstr = unixvmpath(cmdstr, dirmaps);
%%
%ssh call won't allow running a background job with "&" and then quitting
%convert to nohup call instead
isbg = regexp(cmdstr,'[ ]&[ ]*$', 'once');
if(~isempty(isbg) && isbg > 0)
    cmdstr = ['nohup ' cmdstr(1:isbg-1) ' > /dev/null 2>&1 &'];
end

%%
status = [];
status_unmapped = [];
try

    if(isempty(vmconn) || ~test_vmconn(vmconn))
        fprintf('Opening ssh connection to VirtualBox...\n');
        vmconn = init_vmconn(connparams);

    end

    if(do_echo)
        fprintf('%s\n',cmdstr);
    end
    [vmconn,status] = ssh2_command(vmconn,cmdstr,false);
    result = vmconn.command_result;

catch err
    %errors might be:
    %1. virtualbox wasn't running: soln = tell them to start it
    %2. ssh connection was reset: soln = restart connection
    %3. some other ssh problem?
    warning(errorstring(err))
    error('Could not connect to VirtualBox. Either VirtualBox is not running or some other ssh problem occurred.');
end

%how can we get this for real?
exitcode = 0;

if(nargout == 1)
    varargout = {status};
    return;
end

%make a copy of result with directories swapped back to original
status_unmapped = status;
for i = 1:numel(status)
    status_unmapped{i} = unixvmpath(status{i},dirmaps,true);
end
  
if(nargout == 2)
    varargout = {status,status_unmapped};
elseif(nargout == 3)
    varargout = {status,status_unmapped,exitcode};
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isconn = test_vmconn(conn)
result = [];
try
    [conn, status] = ssh2_command(conn,'uname',false);
    result = conn.command_result;
catch err
end

isconn = ~isempty(result);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conn = init_vmconn(connparams)

isconn = false;
try
    conn=ssh2_config_publickey(connparams{:});
    isconn = test_vmconn(conn);
catch err
end

if(~isconn)
    conn = [];
end