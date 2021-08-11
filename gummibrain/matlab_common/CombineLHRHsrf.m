function lrnames = CombineLHRHsrf(filenames)

lstring = '_LH_';
rstring = '_RH_';
lrstring = '_LHRH_';

lhf = regexp(filenames,lstring);
numlhf = sum(~cellfun(@isempty,lhf));

c = 0;

lrnames = {};
for f = 1:numel(lhf)
    if(isempty(lhf{f}))
        continue;
    end
    c = c+1;
    lname = filenames{f};
    rname = [lname(1:lhf{f}-1) rstring lname(lhf{f}+numel(lstring):end)];
    lrname = [lname(1:lhf{f}-1) lrstring lname(lhf{f}+numel(lstring):end)];
    

    if(~(exist(lname,'file') && exist(rname,'file')))
        continue;
    end
    
    fprintf('CombineLHRHsrf %d of %d:\n%s\n%s\n\n',c,numlhf,lrname);        

    lrnames{end+1} = lrname;
    lsrf = xff(lname);
    rsrf = xff(rname);
    lrsrf = lsrf.Combine(rsrf);
    lrsrf.SaveAs(lrname);
end

