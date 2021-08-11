% set defaults
set(0,'DefaultFigureInvertHardCopy','off');
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultFigureColor',[1 1 1]);
set(0,'DefaultFigureColormap',gray(64));
set(0,'DefaultLineLineWidth',1);
set(0,'DefaultLineMarkerSize',9);
set(0,'DefaultTextFontSize',10);
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultTextFontName','Helvetica');
set(0,'DefaultAxesFontName','Helvetica');
%set(0,'DefaultFigureToolbar','none');
%set(0,'DefaultFigureMenuBar','none');
fprintf(1,'default font size: 10\n');

% set rand state
rand('state',sum(100*clock));
randn('state',sum(100*clock));

% set format
format long g;

% define path (top is highest priority)!
pth = '';
pth = [pth genpath('/Users/kjamison/Source/knkutils')];
pth = [pth genpath('/Applications/Psychtoolbox')];

% clean up path
pth = regexp(pth,'(.+?):','match');
bad = {'.svn' '.git' 'DNBdata' 'DNBresults'};
isbad = zeros(1,length(pth));
for p=1:length(pth)
  for q=1:length(bad)
    isbad(p) = isbad(p) | ~isempty(strfind(pth{p},bad{q}));
  end
end

% add to path
pth = pth(~isbad);
pth = cat(2,pth{:});
addpath(pth);

% clean up
clear pth bad isbad p q;

% get to the right directory

