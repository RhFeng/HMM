function pst=p2stat(p);

%
% Stationary distribution from the 
% transition matrix
%

d=size(p,1);
ptr=p';
A=zeros(d,d);
A(1:d-1,1:d)=ptr(1:d-1,:)-[eye(d-1) zeros(d-1,1)];
A(d,1:d)=ones(1,d);
bvec=zeros(d,1);
bvec(d,1)=1;

pst=inv(A)*bvec;
