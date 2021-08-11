function outfile = convert_ant(eepfile)

[datadir filename ext] = fileparts(eepfile);
EEG = pop_loadeep(eepfile , 'triggerfile', 'on');
pop_saveset( EEG, 'filename',[filename '.set'],'filepath',datadir);
outfile = [datadir '/' filename '.set'];