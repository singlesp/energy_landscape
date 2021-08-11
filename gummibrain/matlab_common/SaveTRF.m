function SaveTRF(trf,filename)

fid = fopen(filename,'w');
fprintf(fid,'\nFileVersion: %d\n\nDataFormat: %s\n\n',trf.FileVersion,trf.DataFormat);
fprintf(fid,'%20.16f%20.16f%20.16f%20.16f\n',trf.TFMatrix);
fprintf(fid,'\nTransformationType: %d\nCoordinateSystem: %d\n\n',trf.TransformationType,trf.CoordinateSystem)
fprintf(fid,'SourceFile: "%s"\nTargetFile: "%s"\n',trf.SourceFile,trf.TargetFile)
fclose(fid);
