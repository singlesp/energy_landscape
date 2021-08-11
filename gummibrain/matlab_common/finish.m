try
    savehistory
catch
end

global vmconn;
if(ispc && ~isempty(vmconn))
    try
        ssh2_close(vmconn);
    catch err
    end
end
