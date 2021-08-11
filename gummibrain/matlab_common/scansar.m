%clc;

topdir = '/export/raid2/dicom/20120601-ST001-Jamison2'
%topdir = '/export/raid2/dicom/20120525-ST001-Roy'
%topdir = '/export/raid2/dicom/20120526-ST001-Jamison'
%topdir = '/export/raid2/dicom/20120602-ST001-Jamison6'

scandirs = dir(topdir);
for d = 1:numel(scandirs)
	if(~scandirs(d).isdir)
		continue;
	end
	dcmfiles = dir([topdir '/' scandirs(d).name '/*.dcm']);
	if(isempty(dcmfiles))
		continue;
	end
	dcmfile = [topdir '/' scandirs(d).name '/' dcmfiles(1).name];
	try
		dcminfo = dicominfo(dcmfile);
		if(isfield(dcminfo,'ImageComments'))
			cmt = dcminfo.ImageComments;
		else
			cmt = '';
		end
		
		%sestr = regexp(dcmfiles(1).name,'^MR-ST[0-9]+-SE([0-9]+)-[0-9]','tokens');
		sestr = regexp(scandirs(d).name,'^MR-SE([0-9]+)-(.*$)','tokens');
		seriesnum = sscanf(sestr{1}{1},'%d');
		scanname = sestr{1}{2};
		
		scantype = '';
		switch scanname
      case 'localizer'
        continue;
      case 'cmrr_mbep2d_bold_3x3x3'
        scantype = 'EPI';
      case 'eja_ep2d_bold_2x2x2__MB4'
        if(isfield(dcminfo,'ImageComments'))
          mbtok = regexp(dcminfo.ImageComments,'(MB[0-9]+)','tokens');
          scantype = mbtok{1}{1};
         else
          scantype = 'MB1';
         end
      case 'eja_ep2d_bold_2x2x2__MB4_SBRef'
        continue;
      otherwise
        scantype = scanname;
		end
		fprintf('%04d %-10s %-1.4f %3d %5d %5d %-7.1f %s\n',seriesnum, scantype, dcminfo.SAR, dcminfo.FlipAngle, dcminfo.RepetitionTime, numel(dcmfiles), dcminfo.PatientWeight, scandirs(d).name);
		%fprintf('SAR=%-1.4f TR=%-6d N=%-5d kg=%-4.1f %-24s %s\n',dcminfo.SAR, dcminfo.RepetitionTime, numel(dcmfiles), dcminfo.PatientWeight, %cmt, scandirs(d).name);
		
	catch
		warning('error dicominfo: %s', dcmfile);
		lasterror
	end
end