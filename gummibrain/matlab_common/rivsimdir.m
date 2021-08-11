function d = rivsimdir

if(isunix)
    d = '/home/kjamison/rivalry_eegmri_prism/data';
    %d = '/home/jeet/rivalry_eegmri_prism/data'
elseif(ispc)
    compstr = computerid;
    if(strcmpi(compstr,'workpc64new'))
        d = 'd:/data/rivalry_eegmri_prism';
    else
        d = 'e:/research/rivalry_eegmri_prism';
    end
end
