function hfig = PreviewGUI(varargin)
set(0,'ShowHiddenHandles','off');

zoomfactor = 1;

label_bgcolor = [1 .7 .7];
label_fgcolor = [0 0 0];

if(~iscell(varargin{1}) && ishandle(varargin{1}))
    hfig = varargin{1};
    filenames = varargin{2};
    varargin = varargin(3:end);
else
    hfig = [];
    filenames = varargin{1};
    varargin = varargin(2:end);
end

if(~iscell(filenames) && ischar(filenames))
    filenames = {filenames};
end
if(~iscell(filenames))
    return;
end

maintain_scrollbar_position = false;
if(~isempty(varargin) && islogical(varargin{1}))
    maintain_scrollbar_position = varargin{1};
end

filepref = common_prefix(justfilename(filenames));
filesuff = common_suffix(justfilename(filenames,true));
if(numel(filenames) == 1)
    filesuff = '';
end

imgfiletypes = '(png|jpg|jpeg|gif|bmp)$';
txtfiletypes = '(txt|log|m)$';
matfiletypes = '(mat)$';
zipfiletypes = '(gz)$';

label_fontsize = 9;

%bgcolor = [1 1 1];
bgcolor = get(0,'DefaultFigureColor');
txtfile_color = [1 1 1];
matfile_color = .8*[1 1 1];
    
if(isempty(hfig))
    sufftxt = '';
    if(~isempty(filesuff))
        sufftxt = sprintf('...%s',filesuff);
    end
    hfig = figure('color',bgcolor,'NumberTitle','off','Name',sprintf('Preview: %s%s',filepref,sufftxt));
    maintain_scrollbar_position = false;
end



slider_position = 1;
if(maintain_scrollbar_position)
    hslide = findobj(hfig,'style','slider','tag','PreviewSlider');
    if(~isempty(hslide))
        slider_position = get(hslide,'value');
    end
end
%current_children = get(hfig,'children');
%delete(current_children);
clearchildren(hfig);


outerpanel = uipanel('parent',hfig,'backgroundcolor',bgcolor,'tag','previewpanel');
lblWait = uicontrol('parent',hfig,'style','text','string','Reading files...',...
    'horizontalalignment','center','units','normalized','position',[0 0 1 .5]);

set(outerpanel,'visible','off');
innerpanel = uipanel('parent',outerpanel,'backgroundcolor',bgcolor,'tag','innerpanel');
set(outerpanel,'position',[0 .05 .95 .875],'bordertype','none');
set(innerpanel,'units','pixels','position',[0 0 100 100],'bordertype','none');

zoompanel = uipanel('parent',hfig,'units','normalized','position',[0 .925 1 .075],'bordertype','none');    
lblZoom = uicontrol('parent',zoompanel,'style','text','units','normalized','position',[0 .5 1 .5],...
    'String','zoom','hittest','off','enable','inactive');
sliderZoom = uicontrol('parent',zoompanel,'style','slider','units','normalized','position',[0 0 1 .5],...
    'callback',@sliderZoom_callback,'min',.5,'max',4,'value',1,'sliderstep',[.1 .5]);
set(lblZoom,'String',sprintf('zoom: %.0f%%',round(zoomfactor*100)));

hfig_parent = GetParentFigure(hfig);
parent_pos = get(hfig_parent,'position');
hu = get(hfig,'units');
set(hfig,'units','pixels');
hfig_pos = get(hfig,'position');
set(hfig,'units',hu);
guiratio = parent_pos./hfig_pos;

panelVert = uipanel('parent',hfig,'bordertype','none',...
    'units','normalized','position',[.95 0 .05 1]);
sliderVert = uicontrol('style','slider','parent',panelVert,...
    'units','normalized','position',[0 0 1 1],'value',slider_position,'tag','PreviewSlider');


panelHorz = uipanel('parent',hfig,'bordertype','none',...
    'units','normalized','position',[0 0 .95 .05]);
sliderHorz = uicontrol('style','slider','parent',panelHorz,...
    'units','normalized','position',[0 0 1 1],'value',0);

set(sliderVert,'callback',{@slider_update,sliderVert,sliderHorz,innerpanel});
set(sliderHorz,'callback',{@slider_update,sliderVert,sliderHorz,innerpanel});

%set(panelHorz,'ResizeFcn',{@preview_resize,panelVert,panelHorz});

jvscroll = findjobj(sliderVert);
jvscroll.MouseDraggedCallback = {@slider_update,sliderVert,sliderHorz,innerpanel};
jhscroll = findjobj(sliderHorz);
jhscroll.MouseDraggedCallback = {@slider_update,sliderVert,sliderHorz,innerpanel};

outersize = getpixelposition(outerpanel);

% 
% resizefcn = get(hfig,'ResizeFcn');
% if(iscell(resizefcn))
%     resizefcn{end+1} = {@preview_resize sliderVert sliderHorz} ;
% else
%     resizefcn = {@preview_resize sliderVert sliderHorz};
% end
% set(hfig,'ResizeFcn',resizefcn);

hpadding = 10;
vpadding = 10;
textbox_width = outersize(3)-2*hpadding;
if(numel(filenames) == 1)
    textbox_height = outersize(4) - 2*vpadding - 4*label_fontsize;
else
    textbox_height = outersize(4)/2 - 2*vpadding;
end
textbox_height = textbox_height * zoomfactor;

current_bottom = vpadding;
current_right = 0;
newax = [];
for i = numel(filenames):-1:1
    if(i == 1 || ~mod(i,5))
        set(lblWait,'string',sprintf('%.0f%%...',100*(numel(filenames)-i+1)/numel(filenames)));
    end
    [filedir,~,ext] = fileparts(filenames{i});
    file_to_load = filenames{i};
    
    if(~isempty(regexpi(ext,zipfiletypes)))
        gzip_output = gunzip(filenames{i},tempdir);
        file_to_load = gzip_output{1};
        [d,~,ext] = fileparts(file_to_load);
    end
    
    cmd_for_file = ['xdg-open ' filedir];
    
    if(~isempty(regexpi(ext,imgfiletypes)))
        cmd_for_file = ['xdg-open ' file_to_load];
        img = imread(file_to_load);
        
        newax(end+1) = axes('parent',innerpanel,'units','pixels',...
            'position',[hpadding current_bottom size(img,2)*zoomfactor size(img,1)*zoomfactor]);
        imshow(img,'Parent',newax(end));

    elseif(~isempty(regexpi(ext,matfiletypes)))
        M = load(file_to_load);
        textlines = printstruct(M);
        newax(end+1) = uicontrol('parent',innerpanel,'style','edit','units','pixels',...
            'position',[hpadding current_bottom textbox_width textbox_height],...
            'horizontalalignment','left','max',100,'backgroundcolor',matfile_color,'string',textlines);
        
    elseif(~isempty(regexpi(ext,txtfiletypes)))
        cmd_for_file = ['gedit ' file_to_load];
        
        fid = fopen(file_to_load,'r');
        textlines = {};
        str = fgetl(fid);
        while ischar(str)
            textlines{end+1} = str;
            str = fgetl(fid);
        end
        fclose(fid);
        
        newax(end+1) = uicontrol('parent',innerpanel,'style','edit','units','pixels',...
            'position',[hpadding current_bottom textbox_width textbox_height],...
            'horizontalalignment','left','max',100,'backgroundcolor',txtfile_color,'string',textlines);

    else
        warning('unknown file type: %s', filenames{i});
    end
    
    current_right = max(current_right, element(get(newax(end),'position'),3) + hpadding);
    current_bottom = current_bottom + element(get(newax(end),'position'),4);
    
    textheight = 30;
    liststr = sprintf('(%d/%d) %s',i,numel(filenames),justfilename(filenames{i}));
    pos = [0 current_bottom textbox_width textheight];
    newax(end+1) = uicontrol('parent',innerpanel,'style','text','string',liststr,...
        'fontsize',label_fontsize,'horizontalalignment','left','position',pos,'backgroundcolor',label_bgcolor,...
        'buttondownfcn',@label_click,...
        'userdata',struct('dir',filedir,'filename',filenames{i},...
        'loaded_file',file_to_load,'cmd',cmd_for_file));
    [liststr_wrapped newpos] = textwrap_nospaces(newax(end),{liststr});
    set(newax(end),'string',liststr_wrapped,'position',newpos);
    set(newax,'tag','preview');

    current_bottom = current_bottom + vpadding + element(get(newax(end),'position'),4);
end

for i = 1:numel(newax)
    p = get(newax(i),'position');
    u = get(newax(i),'userdata');
    if(isempty(u))
        u = struct('origposition',p);
    else
        u = mergestruct(u,'origposition',p);
    end
    set(newax(i),'userdata',u);
end
    
set(sliderZoom,'userdata',fillstruct(innerpanel));
set(innerpanel,'position',[0 0 current_right current_bottom]);

set(innerpanel,'userdata',fillstruct(innerpanel,outerpanel,panelVert,...
    panelHorz,sliderVert,sliderHorz,vpadding,hpadding,label_fontsize,zoompanel,lblZoom));

delete(lblWait);
set(outerpanel,'visible','on');

%slider_update([],[],sliderVert,sliderHorz,innerpanel);

set(outerpanel,'ResizeFcn',{@preview_resize,panelVert,panelHorz},'interruptible','off');
preview_resize(outerpanel,[],panelVert,panelHorz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_update(src,eventdata,hvert,hhorz,hpanel)
vval = get(hvert,'value');
hval = get(hhorz,'value');
parent = get(hpanel,'parent');


%set(parent,'units','pixels');
parentsize = getpixelposition(parent);
set(parent,'units','normalized');
mysize = get(hpanel,'position');

hpos = mysize(1);
vpos = mysize(2);

if(parentsize(4) < mysize(4))
    vpos = vval*(parentsize(4)-mysize(4));
else
    vpos = 0;
end

if(parentsize(3) < mysize(3))
    hpos = hval*(parentsize(3)-mysize(3));
end

set(hpanel,'position',[hpos vpos mysize(3) mysize(4)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_click(src,eventdata)
M = get(src,'userdata');
if(isstruct(M))
    system(M.cmd);
end
%open(M.loaded_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderZoom_callback(hObj,event)

zoomfactor = get(hObj,'value');
hpanel = mainpanel(hObj);
ch= findobj(hpanel,'tag','preview');

A = get(hpanel,'userdata');

panelpos = get(hpanel,'position');
parpos = getpixelposition(get(hpanel,'parent'));
panelpos(3) = parpos(3);
set(hpanel,'position',parpos);
vpadding = A.vpadding;
curpos = 0;
curright = 0;
for i = numel(ch):-1:1
    u = get(ch(i),'userdata');
    pos = u.origposition;
    pos(3:4) = pos(3:4)*zoomfactor;
    pos(2) = curpos;
    set(ch(i),'position',pos);
    if(strcmpi(get(ch(i),'type'),'uicontrol') && strcmpi(get(ch(i),'style'),'text'))
        %set(ch(i),'fontsize',A.label_fontsize*zoomfactor);
        pos(3) = panelpos(3);
        set(ch(i),'position',pos);

        str = get(ch(i),'string');
        [str_wrapped pos] = textwrap_nospaces(ch(i),{strcat(str{:})});
        set(ch(i),'string',str_wrapped,'position',pos);
        
    end
    curright = max(curright,pos(1)+pos(3));
    curpos = curpos + vpadding + pos(4);
end
sz = get(hpanel,'position');
sz(3) = max(curright,panelpos(3));
sz(4) = curpos + vpadding;
set(hpanel,'position',sz);
set(A.lblZoom,'String',sprintf('zoom: %.0f%%',round(zoomfactor*100)));
%preview_resize(A.innerpanel,[],A.panelVert,A.panelHorz);
%slider_update([],[],A.sliderVert,A.sliderHorz,A.innerpanel);
%set(hpanel
%prevobj = findobj(get(hObj,'parent'),'tag','preview')
%get(hObj,'value')
%get(prevobj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preview_resize(src,eventdata, panelVert, panelHorz)

sliderpx = 20;


hparent = get(src,'parent');
innerpanel = mainpanel(src);
outerpanel = get(innerpanel,'parent');

A = get(innerpanel,'userdata');

parentpos = getpixelposition(hparent);


vpos = [parentpos(3)-sliderpx sliderpx sliderpx parentpos(4)-sliderpx];
hpos = [0 0 parentpos(3)-sliderpx sliderpx];

innerpos = [0 sliderpx parentpos(3)-sliderpx parentpos(4)];

zpos = getpixelposition(A.zoompanel);
zpos(3) = parentpos(3)-sliderpx;
setpixelposition(A.zoompanel,zpos);

%get(outerpanel)
outerpos = [0 sliderpx parentpos(3)-sliderpx parentpos(4)-zpos(4)-sliderpx];
%setpixelposition(outerpanel,outerpos);

setpixelposition(panelVert,vpos);
setpixelposition(panelHorz,hpos);
if(innerpanel > 0)
    setpixelposition(innerpanel,innerpos);
end

txtbox = findobj(hparent,'type','uicontrol','style','edit');
for t = 1:numel(txtbox)
    txtpos = getpixelposition(txtbox(t));
    txtpos(3) = parentpos(3)-sliderpx-2*txtpos(1);
    setpixelposition(txtbox(t),txtpos);
end

lbl = findobj(innerpanel,'tag','preview');
if(~isempty(lbl))
    px1 = cellfun(@getpixelposition,num2cell(lbl),'uniformoutput',false);
    mx = max(cellfun(@(x)(x(1)+x(3)),px1));
    my = max(cellfun(@(x)(x(2)+x(4)),px1));
    
    px2 = getpixelposition(innerpanel);
    if(px2(3) < mx)
        px2(3) = mx;
        
    end
%     hcontainer = get(lbl(1),'parent');
%     hconparent = get(hcontainer,'parent');
%     
%     px1 = getpixelposition(hcontainer);
%     px2 = getpixelposition(hconparent);
%     
    if(px2(4) < my)
        px2(4) = my;
    end
    setpixelposition(innerpanel,px2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hpanel = mainpanel(hObj)
curtag = get(hObj,'tag');
while(~strcmpi(curtag,'previewpanel'))
    hObj = get(hObj,'parent');
    if(hObj == 0)
        hObj = gcbf;
        break;
    end
    curtag = get(hObj,'tag');
end

hpanel = findobj(hObj,'tag','innerpanel');
    
