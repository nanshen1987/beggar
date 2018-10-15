
% num = 50; % num¸öËæ»úµã
% Rand1 = randi([-1,1],num,3); %ÔëÉù·¶Î§
% 
% X=(1:0.5:0.5*(num+1))'+Rand1(:,1);
% Y=(1:0.5:0.5*(num+1))'+Rand1(:,2);
% Z=(1:0.5:0.5*(num+1))'+Rand1(:,3);

close all
X=gpsData.X(50:90);
Y=gpsData.Y(50:90);
Z=gpsData.Z(50:90);
x=X'-mean(X);
y=Y'-mean(Y);
z=Z'-mean(Z);




F=[z;1 1 1 1 1 1 1 1 1 1 1];
M=F*F';
N=F*x';
O=F*y';
A=(M\N)';
B=(M\O)';
x1=A(1)*z+A(2);
y1=B(1)*z+B(2);
z1=z;
plot3(x1,y1,z1,'b',x,y,z,'o');


% [a,resnorm]=fit_line([1 1 1 1 1 1],[X,Y,Z])
% plot3(X,Y,Z,'r+')
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold on
% t=-30:30;
% plot3(a(1)+t*a(4),a(2)+t*a(5),a(3)+t*a(6));
