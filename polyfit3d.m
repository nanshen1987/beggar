function [p,error] = polyfit3d(X,Y,Z)
%POLYFIT3D 此处显示有关此函数的摘要
%   此处显示详细说明
[n,~]=size(X);
B=zeros(n*2,4);
L=zeros(n*2,1);
for i=1:n
    B(i,1)=Z(i);
    B(i,2)=1;
    B(i+1,3)=Z(i);
    B(i+1,4)=1;
    L(i)=X(i);
    L(i+1)=Y(i);
end
error=1;
p=(B'*B)\(B'*L);
figure 

procX=p(1)*Z+p(2);
procY=p(3)*Z+p(4);


plot3(procX,procY,Z,'r+',X,Y,Z,'b-');
end

