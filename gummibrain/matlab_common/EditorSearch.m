function EditorSearch

global EditorSearchData;

if(isempty(EditorSearchData))
    do_comments = false;
    do_word = false;
    do_wrap = true;
    do_case = false;
    searchtxt = '';
    replacetxt = '';
else
    dumpstruct(EditorSearchData);
end

%clh;
hfig = findall(0,'tag','figEditorSearch');
if(isempty(hfig))
    hfig = hgload('EditorSearch_gui.fig',struct('visible','off'));
end
set(hfig,'visible','off');
uiobj = FigureTagChildren(hfig);
dumpstruct(uiobj);

searchin_fmt = {'Editor - Current Selection',...
    'Editor - Current File (%s)',...
    'Editor - All Files'};
searchin_val = {'cursel','curfile','allfile'};

file = matlab.desktop.editor.getActive;
curfile = file.Filename;
seltxt = file.SelectedText;
if(~isempty(seltxt))
    seltxt = regexp(seltxt,'[\r\n]+','split');
    seltxt = seltxt{1};
end

if(isempty(seltxt))
    set(uiobj.popupSearch,'value',find(strcmp(searchin_val,'curfile')));
else
    searchtxt = seltxt;
    set(uiobj.popupSearch,'value',find(strcmp(searchin_val,'cursel')));
end
if(~isempty(searchtxt))
    set(uiobj.txtFind,'string',searchtxt);
end
uicontrol(uiobj.txtFind);
if(~isempty(replacetxt))
    set(uiobj.txtReplace,'string',replacetxt);
end
G = fillstruct(uiobj,searchin_fmt,searchin_val,searchtxt,replacetxt,...
    do_comments,do_wrap,do_word,do_case,curfile);

setappdata(hfig,'appdata',G);
update_gui(hfig);

fn = fieldnames(uiobj);
for i = 1:numel(fn)
    obj = uiobj.(fn{i});
    if(isfield(get(obj),'Callback'))
        set(obj,'callback',{@editsearch_callback,fn{i}});
    end
end
set(hfig,'units','pixels','CloseRequestFcn',{@editsearch_callback,'closefig'});
p_default = get(0,'defaultfigureposition');
p = get(hfig,'position');
p(1:2) = p_default(1:2);
set(hfig,'units','pixels','position',p,'visible','on');
EditorSearchData = G;
WinOnTop(hfig,true);
%%

function update_gui(fig)
G = getappdata(fig,'appdata');
dumpstruct(G);
for i = 1:numel(searchin_fmt)
    if(any(searchin_fmt{i}=='%'))
        searchin_fmt{i} = sprintf(searchin_fmt{i},justfilename(curfile));
    end
end
set(uiobj.popupSearch,'String',searchin_fmt);

set(uiobj.chkComments,'value',bool2check(do_comments));
set(uiobj.chkWrap,'value',bool2check(do_wrap));
set(uiobj.chkWord,'value',bool2check(do_word));
set(uiobj.chkCase,'value',bool2check(do_case));
%%
function c = bool2check(b)
if(b)
    %c = 'on';
    c = 1;
else
    %c = 'off';
    c = 0;
end
%%
function b = check2bool(c)
if(ishandle(c))
    %c = get(c,'checked');
    c = get(c,'value');
end
%b = strcmpi(c,'on');
b = c > 0;

%%
function r = editsearch_callback(gcbo,ev,param)
global EditorSearchData;

fig = gcf;
G = getappdata(fig,'appdata');
if(isfield(G,'mainfig'))
    fig = G.mainfig;
    G = getappdata(fig,'appdata');
end
dumpstruct(G);
dumpstruct(uiobj);


paramvalues = {};
if(iscell(param))
    paramvalues = param(2:end);
    param = param{1};
end

searchtxt= get(uiobj.txtFind,'string');
replacetxt = get(uiobj.txtReplace,'string');
searchtype = searchin_val{get(uiobj.popupSearch,'value')};

switch param

    case 'cmdClose'

        G.searchtxt = searchtxt;
        G.replacetxt = replacetxt;
        EditorSearchData = G;
        close(fig);
        return;
    case 'closefig'

        G.searchtxt = searchtxt;
        G.replacetxt = replacetxt;
        EditorSearchData = G;
        delete(fig);
        return;
    case 'popupSearch'
        
    case 'cmdFindPrev'
        searchdir = -1;
        replacestr = '';
        r = perform_search(searchtype, searchdir, searchtxt, replacetxt, false, do_wrap, do_word, do_case, do_comments);
    case 'cmdFind'
        searchdir = 1;
        replacestr = '';
        r = perform_search(searchtype, searchdir, searchtxt, replacetxt, false, do_wrap, do_word, do_case, do_comments);
    case 'cmdReplace'
        searchdir = 1;
        r = perform_search(searchtype, searchdir, searchtxt, replacetxt, true, do_wrap, do_word, do_case, do_comments);
    case 'cmdReplaceAll'
        searchdir = 0;
        r = perform_search(searchtype, searchdir, searchtxt, replacetxt, true, do_wrap, do_word, do_case, do_comments);        
    case 'chkComments'
        G.do_comments = check2bool(uiobj.(param));
    case 'chkWord'
        G.do_word = check2bool(uiobj.(param));
    case 'chkWrap'
        G.do_wrap = check2bool(uiobj.(param));
    case 'chkCase'
        G.do_case = check2bool(uiobj.(param));
        
end
setappdata(fig,'appdata',G);

%%
function r = perform_search(searchin_val, searchdir, searchstr, replacestr, do_replace, do_wrap, do_word, do_case, do_comments)
r = 1;

if(isempty(searchstr))
    return;
end

do_sel = false;
switch searchin_val
    case 'cursel'
        file = matlab.desktop.editor.getActive;
        do_sel = true;
        if(searchdir < 0)
            searchdir = 1;
            warning('previous search in selection not supported');
        end
    case 'curfile'
        file = matlab.desktop.editor.getActive;
    case 'allfile'
        file = matlab.desktop.editor.getAll;
end

if(iscell(searchstr))
    if(numel(searchstr) == 1)
        searchstr = searchstr{1};
    else
        searchstr = sprintf('%s\n',searchstr{:});
        searchstr = searchstr(1:end-1);
    end
end

searchstr_orig = searchstr;
wordchars = '[^A-Za-z0-9]';
if(do_word)
    searchstr = [wordchars searchstr wordchars];
end

if(~do_case)
    regfcn = @regexpi;
    casestr = {'ignorecase'};
else
    regfcn = @regexp;
    casestr = {};
end
    
file_restarted = false;

while(~isempty(file))
    if(isempty(file(1).Text))
        continue;
    end
    
    txt = file(1).Text;
    if(~file_restarted)
        selline = sort(file(1).Selection([1 3]));
        selchar = sort(file(1).Selection([2 4]));
        sel1 = matlab.desktop.editor.positionInLineToIndex(file(1),selline(1),selchar(1));
        sel2 = matlab.desktop.editor.positionInLineToIndex(file(1),selline(2),selchar(2));
        if(searchdir > 0)
            searchrange = [sel1+1 numel(txt)];
        elseif(searchdir < 0)
            searchrange = [1 sel1];
        elseif(searchdir == 0)
            searchrange = [1 numel(txt)];
        end
    end
    file_restarted = false;

    
    if(do_sel && ~isempty(file(1).SelectedText))
        if(sel2 > sel1)
            searchrange = [sel1+1 sel2];
        end
    end

    ridx = regfcn(txt(searchrange(1):searchrange(2)),searchstr);
    r_found = [];
    
    if(~isempty(ridx))
        if(do_word)
            ridx = ridx + 1;
        end
        txtlines =  matlab.desktop.editor.textToLines(txt);
        if(searchdir < 0)
            rrange = numel(ridx):-1:1;
        else
            rrange = 1:numel(ridx);
        end
        for i = rrange
            [l1, c1] = matlab.desktop.editor.indexToPositionInLine(file(1),searchrange(1)+ridx(i)-1);
            if(do_comments || isempty(regexp(txtlines{l1},'^\s+%')))
                r_found = [r_found ridx(i)];
                if(searchdir ~= 0)
                    break;
                end
            end
        end
    end
    if(isempty(r_found))
        if(do_wrap && searchdir ~= 0)
            file_restarted = true;
            searchrange = [1 numel(txt)];
            continue;
        else
            file = file(2:end);
            continue;
        end
    end
    if(searchdir == 0)
        r_found = sort(r_found,'descend');
    end
    
    [l1, c1] = matlab.desktop.editor.indexToPositionInLine(file(1),searchrange(1)+r_found(1)-1);

    if(do_replace)
        for i = 1:numel(r_found)
            r_start = searchrange(1)+r_found(i)-1;
            pre_win = 1:r_start-1;
            post_win = (r_start+numel(searchstr_orig)):numel(txt);
            txt = [txt(pre_win) replacestr txt(post_win)];
        end

        file(1).makeActive();
        file(1).Text = txt;
        if(searchdir == 0)
            file(1).Selection = [selline(1) selchar(1) selline(2) selchar(2)];
        else
            file(1).Selection = [l1 c1 l1 c1+numel(replacestr)];
        end
        
    else
        %file(1).goToPositionInLine(l1,c1);
        file(1).makeActive();
        file(1).Selection = [l1 c1 l1 c1+numel(searchstr_orig)];
        %
    end
    file = file(2:end);
end


