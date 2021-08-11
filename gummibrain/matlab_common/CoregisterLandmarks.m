function pos2 = CoregisterLandmarks(pos1, landmarks1, landmarks2)
%pos1 = 3xN position
%landmarks1 = 3x3 column vectors with landmarks for pos1
%landmarks2 = 3x3 column vectors with new landmarks
%return: 3xN position relative to landmarks2

A = landmarks1;
B = landmarks2;

A_perp = cross(A(:,2)-A(:,1), A(:,3)-A(:,1));
A_perp = A_perp/norm(A_perp);

B_perp = cross(B(:,2)-B(:,1),B(:,3)-B(:,1));
B_perp = B_perp/norm(B_perp);

A_size = tril(distance(A,A));
B_size = tril(distance(B,B));
A_size = mean(A_size(:));
B_size = mean(B_size(:));

A = [A A_size*A_perp+A(:,1)];
B = [B B_size*B_perp+B(:,1)];

A = [A; ones(1,size(A,2))];
B = [B; ones(1,size(B,2))];

AtoB = B/A;

pos1 = [pos1; ones(1,size(pos1,2))];

pos2 = AtoB*pos1;
pos2 = pos2(1:3,:);