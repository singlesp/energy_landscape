function d = freesurferdir

if(ispc)
    d = sprintf('%s/freesurfer',researchdir);
else
    [~,d] = system('echo $SUBJECTS_DIR');
end

%d = [d '_xinyang/DESPIKE'];
