function Vfs = mesh_fieldsign(Veccen,Vpolar,pos,Ne,smoothiter,eccen_offset,polar_offset)
if(~exist('eccen_offset','var') || isempty(eccen_offset))
    eccen_offset = 0;
end
if(~exist('polar_offset','var') || isempty(polar_offset))
    polar_offset = 0;
end

%Veccen = mesh_diffuse(Veccen,Ne,smoothiter);
%Vpolar = mesh_diffuse(Vpolar,Ne,smoothiter);

A = mesh_diffuse([real(Veccen) imag(Veccen) real(Vpolar) imag(Vpolar)],...
        Ne,smoothiter);
Veccen = complex(A(:,1),A(:,2));
Vpolar = complex(A(:,3),A(:,4));

Aeccen = angle(Veccen) + eccen_offset;
Apolar = angle(Vpolar) + polar_offset;

Vfa = zeros(size(Veccen));

%spherical coordinates
if(size(pos,2) == 3 && ~all(pos(:,3)==pos(1,3)))
    th = atan2(pos(:,2),pos(:,1));
    d = sqrt(pos(:,1).*pos(:,1)+pos(:,2).*pos(:,2));
    ph = atan2(d,pos(:,3));
    pos = [-th ph];
end

Vfs = mesh_fieldsign_double(Aeccen,Apolar,pos,Ne);
return;
%%
 for i = 1:size(Veccen,1)
     if(isempty(Ne{i}))
         continue;
     end
     %circsubtract
     dR = wrapToPi(Aeccen(i) - Aeccen(Ne{i}));
     dTh = wrapToPi(Apolar(i) - Apolar(Ne{i}));
     
     if(all(isnan(dR)))
         continue;
     end
     %dp = bsxfun(@minus,pos(i,:),pos(Ne{i},:));
     %du = dp(:,1);
     %dv = dp(:,2);
     du = pos(i,1)-pos(Ne{i},1);
     dv = pos(i,2)-pos(Ne{i},2);
     
     y = [dR dTh];
     x = [du dv ones(size(du))];
     
     %Ax = B ----> x = pinv(A)*B .... min(Ax - B)
%			%so A = [du dv 1], B = [dR dTh]
     %y = xB+n
     %B = pinv(x)*y = x\y;
     B = x\y;

     Vfa(i) = det(B(1:2,1:2));
     %%if(i<=3)
     %%    fprintf('%20.3f %7.3f\n',B.');
     %%    fprintf('det=%.3f\n',Vfa(i));
     %%end
 end
%%toc

Vfs = Vfa;

%Vfs = sign(Vfa);
%Vstat = 1000*abs(Vfa);
%Vmask = sqrt(abs(Veccen).*abs(Vpolar)); %geom mean