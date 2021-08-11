function [img, ser, mrprot] = parse_siemens_shadow(dcm)
% [img, ser, mrprot] = parse_siemens_shadow(dcm)
% function to parse siemens numaris 4 shadow data
% returns three structs with image, series header, mrprot info
% does not work with arrayed dcm()

if (size(dcm,2) > 1)
    error('parse_siemens_shadow does not work on arrayed dicominfo data!')
end

ver_string = dcm.Private_0029_1008;
csa_string = dcm.Private_0029_10xx_Creator;

if (strcmp(ver_string,'IMAGE NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA HEADER'))
        img = parse_shadow_func(dcm.Private_0029_1010);
        ser = parse_shadow_func(dcm.Private_0029_1020);
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
elseif (strcmp(ver_string,'SPEC NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA NON-IMAGE'))
        if isfield(dcm,'Private_0029_1210')
            img = parse_shadow_func(dcm.Private_0029_1210);
            ser = parse_shadow_func(dcm.Private_0029_1220);
        else %VB13
            img = parse_shadow_func(dcm.Private_0029_1110);
            ser = parse_shadow_func(dcm.Private_0029_1120);
        end
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
else
    error('shadow: Unknown/invalid NUMARIS version: %s',ver_string);
end

% now parse the mrprotocol
tmp_fn = tempname;
fp = fopen(tmp_fn,'w+');
if isfield(ser, 'MrPhoenixProtocol') % VB13
    MrProtocol = char(ser.MrPhoenixProtocol);
else
    MrProtocol = char(ser.MrProtocol);
end
spos = strfind(MrProtocol,'### ASCCONV BEGIN ###');
epos = strfind(MrProtocol,'### ASCCONV END ###');
MrProtocol = MrProtocol(spos+22:epos-2);
fwrite(fp,MrProtocol,'char');
frewind(fp);
mrprot = parse_mrprot(fp);
fclose(fp);
delete(tmp_fn);

%--------------------------------------------------------------------------

function hdr = parse_shadow_func(dcm)
% internal function to parse shadow header

% dump everything to a temp file
tmp_fn = tempname;

% open it using little endian ordering, so this should work on any machine
% (input data is uint8, so write is always ok)
fp = fopen(tmp_fn,'w+','ieee-le');
fwrite(fp,dcm);
frewind(fp);

%fprintf('\nReading IMAGE header\n');
%img_ver = fread(fp,4,'*char')';                  % version string?  SV10
%int1 = fread(fp,1,'int32');                        % unknown (chars 1,2,3,4)
fseek(fp,8,'cof');
nelem = fread(fp,1,'int32');                       % # of elements
fseek(fp,4,'cof'); %int2 = fread(fp,1,'int32');                        % unknown (77)
%fprintf('Found %d elements\n', nelem);
for y=1:nelem
    %data_start_pos = ftell(fp);
%rrr    tag = c_str(fread(fp,64,'*char')');
    tag = c_str(fread(fp,64,'uint8')');
    
    tag = strrep(tag,'-','_'); % remove invalid chars from field name
    vm = fread(fp,1,'int32');
%rrr    vr = c_str(fread(fp,4,'*char')');
    vr = c_str(fread(fp,4,'uint8')');    
    
    fseek(fp,4,'cof'); %SyngoDT = fread(fp,1,'int32');                 % SyngoDT
    NoOfItems = fread(fp,1,'int32') ;              % NoOfItems

 %   keyboard

    if (vm == 0) % this can happen in spectroscopy files, VB13 image files
        %vm = 6;
        vm = NoOfItems;
    end
    
    %str_data = '';
    
    if (NoOfItems > 1)
        fseek(fp,4,'cof'); %int3 = fread(fp,1,'int32'); % unknown (77)

        for z=1:vm
            int_items = fread(fp,4,'int32'); % unknown (77)
            use_fieldwidth = ceil(int_items(4) / 4) * 4;
%rrr            tmp_data = c_str(fread(fp,use_fieldwidth,'*char')');
            tmp_data = c_str(fread(fp,use_fieldwidth,'uint8')');
            %str_data = [str_data tmp_data];
            %if (z < vm), str_data = [str_data '\']; end

            switch vr
                case {'AE','AS','CS','DA','DT','LO','LT','OB','OW','PN','SH','SQ','ST','TM','UI','UN','UT'}
                    % these are string values
                    %fprintf('String VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
                case {'IS','LO','SL','SS','UL','US'}
                    % these are int/long values
                    %fprintf('%s: Int/Long VM %d, VR %s, data = %s, val = %d\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                case {'DS','FL','FD'} 
                    % these are floating point values
                    %fprintf('%s: Float/double VM %d, VR %s, data = %s, val = %.8f\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                otherwise % just assume string
                    %error('Unknown VR = %s found!\n',vr);
                    %fprintf('Unknown VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
            end
        end
    else
        val_data = [];
    end

    junk_len = 16 * (NoOfItems - vm);
    if (NoOfItems < 1)
        junk_len = 4;
    else
        if ( (junk_len < 16) && (junk_len ~= 0) )
            junk_len = 16;
        end
    end

    %data_end_pos = ftell(fp);
    fseek(fp,junk_len,'cof'); %junk_data = fread(fp,junk_len,'*char');
    %data_final_pos = ftell(fp);

%     fprintf('%2d - ''%s''\tVM %d, VR %s, SyngoDT %d, NoOfItems %d, Data',y-1, tag, vm, vr, SyngoDT, NoOfItems);
%     if (size(str_data))
%         fprintf(' ''%s''', str_data);
%     end
%     fprintf('\n');
%     fprintf('data_len: %d  pad_len: %d  total_len: %d\n',data_end_pos-data_start_pos,data_final_pos-data_end_pos,data_final_pos-data_start_pos);
%     fprintf('VR %s\n',vr);

    hdr.(tag) = val_data;
end

fclose(fp);
delete(tmp_fn);


function s = c_str(in)
s = char(in);
