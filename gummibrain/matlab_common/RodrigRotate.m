function [pnew Rmat] = RodrigRotate(a,theta,p)

A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
R = eye(3) + A*sin(theta) + A*A*(1-cos(theta));
R = expm(A*theta);
%k = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]
%R = eye(3) + k*sin(theta) + (1-cos(theta))*k*k

if(nargin > 2 && ~isempty(p))
    pnew = [R*p']';
else
    pnew = [];
end

if(nargout > 1)
    Rmat = R;
end