function ch = FigureTagChildren(fig)

if(~exist('fig','var'))
    fig = gcf;
end

%% new matlab tags other figure elements, so find these and remove them
fig2 = figure('visible','off');
o_default = findobj(fig2);
t_default = get(o_default,'tag');
if(ischar(t_default))
    t_default = {t_default};
end
tt_default = ~cellfun(@isempty,t_default);
o_default = num2cell(o_default(tt_default));
t_default = t_default(tt_default);
close(fig2);
%%

o = findobj(fig);
t = get(o,'tag');
if(ischar(t))
    t = {t};
end
tt = ~cellfun(@isempty,t);
o = num2cell(o(tt));
t = t(tt);

%%
[~, i_new] = setdiff(t,t_default);
i_new(strcmp(t(i_new),'none')) = [];
t = t(i_new);
o = o(i_new);
%%


c = reshape([t o].',[],1);
ch = struct(c{:});
