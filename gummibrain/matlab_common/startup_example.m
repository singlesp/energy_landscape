% %------------ FreeSurfer -----------------------------%
% fshome = getenv('FREESURFER_HOME');
% fsmatlab = sprintf('%s/matlab',fshome);
% if (exist(fsmatlab) == 7)
%     path(path,fsmatlab);
% end
% clear fshome fsmatlab;
% %-----------------------------------------------------%
% 
% %------------ FreeSurfer FAST ------------------------%
% fsfasthome = getenv('FSFAST_HOME');
% fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);
% if (exist(fsfasttoolbox) == 7)
%     path(path,fsfasttoolbox);
% end
% clear fsfasthome fsfasttoolbox;
% %-----------------------------------------------------%
% 
% % so that we can override text.m to avoid annoying bug
% warning('off','MATLAB:dispatcher:nameConflict');

% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup'), PsychStartup; end;

%com.mathworks.mlwidgets.html.HtmlComponentFactory.setBrowserProperty('JxBrowser.BrowserType','Mozilla15');
