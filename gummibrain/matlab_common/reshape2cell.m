function outcell = reshape2cell(list,cell_template)

outcell = deepcopy(cell_template);
reshape2cell_mex(list,outcell);
