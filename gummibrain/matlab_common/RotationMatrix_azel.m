function T = RotationMatrix_azel(azimuth,elevation,tilt)
T=RotationMatrix(-elevation*pi/180,-tilt*pi/180,-azimuth*pi/180);
% if(~exist('tilt','var') || isempty(tilt))
%     tilt=0;
% end
% 
% viewN=numel(azimuth);
% %view_az=azimuth(:)*pi/180 -pi/2;
% %view_el=elevation(:)*pi/180;
% %view_tilt=tilt(:)*pi/180;
% %view_az=azimuth(:)*pi/180 -pi/2;
% %view_el=elevation(:)*pi/180 + pi/2;
% 
% %view_az=azimuth(:)*pi/180;
% %view_el=elevation(:)*pi/180;
% %view_tilt=tilt(:)*pi/180;
% 
% view_az=azimuth(:)*pi/180 + pi/2;
% view_el=elevation(:)*pi/180;% - pi/2;
% view_tilt=tilt(:)*pi/180;
% 
% %only X is flipped:
% % view_az=azimuth(:)*pi/180 -pi/2;
% % view_el=elevation(:)*pi/180 + pi/2;
% % view_tilt=tilt(:)*pi/180;
% 
% view_az0=azimuth(:)*pi/180;
% 
% [x,y,z]=sph2cart(view_az,view_el,1);
% viewpts=[x y z]
% 
% 
% T=zeros(4,4,viewN);
% 
% for i = 1:viewN
% 
%     az_tilt=0;
%     
%     v3=viewpts(i,:);
%     v3=v3/norm(v3)
%     v3(1:2)*v3(1:2)'
%     if(sqrt(v3(1:2)*v3(1:2)') < 1e-3)
%         %basically along the Z-axis
%         v2=cross(v3,[0 1 0]);
%         v2=v2/norm(v2);
%         v1=-cross(v3,v2);
%         v1=v1/norm(v1);
%         
%         if(v3(3)<0)
%             az_tilt=-view_az0; %%%%!!!! ugly...
%             if(abs(az_tilt)<1e-3)
%                 az_tilt=az_tilt+pi;
%             end
%         else
%             az_tilt=view_az0;
%         end
%         
%         
% 
%     else
%         v2=cross(v3,[0 0 1]);
%         v2=v2/norm(v2);
%         v1=cross(v2,v3);
%         v1=v1/norm(v1);
%     end
%     
%     Mtilt=eye(4);
%     Mtilt(1:2,1:2)=[cos(view_tilt(i)+az_tilt) sin(view_tilt(i)+az_tilt); 
%         -sin(view_tilt(i)+az_tilt) cos(view_tilt(i)+az_tilt)];
% 
%     M=eye(4);
%     M(1:3,1:3)=inv([v1(:) v2(:) v3(:)]);
%     
%     T(:,:,i)=Mtilt*M;
% end