function [pass,bias,distance] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins)
%% check condition for allignment
%   nodeSol gpsNode for check
%   imuData ins data
%   gpsWinsize allgin window
%   insGpsRatio ratio  of gps and ins
%   epoch4ins current epoch of ins

limitBias=6;
limitDistance=50;
limitR=2.8641;

[~,nodeSolSize]=size(nodeSol.X);
[imuDataSize,~]=size(imuData.gpsSecond);
%%linear check

XYZ=zeros(gpsWinsize,3);
index=1;
for i=nodeSolSize:-1:(nodeSolSize-gpsWinsize+1)
    XYZ(index,1)=nodeSol.X(i);
    XYZ(index,2)=nodeSol.Y(i);
    XYZ(index,3)=nodeSol.Z(i);
    index=index+1;
end
avgXYZ = mean(XYZ, 1);
subtractedXYZ = bsxfun(@minus, XYZ, avgXYZ);
[~, ~, V] = svd(subtractedXYZ);
direction = V(:, 1);
p0 = avgXYZ;
d = direction;
sumBias=0;
for i=1:gpsWinsize
   biasD=norm(cross(XYZ(i,:)-p0,d))/norm(d);
   sumBias=sumBias+biasD*biasD;
end
bias=sqrt(sumBias/(gpsWinsize-1));
%% distance check
distance=norm(XYZ(1)-XYZ(gpsWinsize));

%% distance bias ratio
dbRatio=distance/bias;


%% decide pass or not
if bias<limitBias && distance>limitDistance 
    pass=1;
else
    pass=0;
end
%% cal b2e
% cal ins start and end index
gpsStartSec=nodeSol.gpsSecond(nodeSolSize-gpsWinsize+1);
gpsEndSec=nodeSol.gpsSecond(nodeSolSize);
insStartIndex=epoch4ins;
insEndIndex=epoch4ins;
if pass
    for i=epoch4ins:imuDataSize
        if imuData.gpsSecond(i)>=gpsEndSec
            insEndIndex=i;
            break;
        end
    end
    for i=epoch4ins:imuDataSize
        if imuData.gpsSecond(i)>=gpsStartSec
            insStartIndex=i;
            break;
        end
    end
%cbe估计值
B4Cbe=zeros((gpsWinsize-1)*3,9);
L4Cbe=zeros((gpsWinsize-1)*3,1);
iter4Cbe=0;
for i=nodeSolSize-gpsWinsize+1:nodeSolSize-1
    L4Cbe(iter4Cbe*3+1)=nodeSol.Vx(i+1)-nodeSol.Vx(i);
    L4Cbe(iter4Cbe*3+2)=nodeSol.Vy(i+1)-nodeSol.Vy(i);
    L4Cbe(iter4Cbe*3+3)=nodeSol.Vy(i+1)-nodeSol.Vz(i);
    sectSt=nodeSol.gpsSecond(i);
    sectEnd=nodeSol.gpsSecond(i+1);
    sumVx=0;
    sumVy=0;
    sumVz=0;
    for j=insStartIndex:insEndIndex-1
        imuStart=imuData.gpsSecond(j);
        if imuStart>=sectSt&&imuStart<sectEnd
            deltaTcbe=imuData.gpsSecond(j+1)-imuStart;
            sumVx=sumVx+(imuData.accX(j)+imuData.accX(j+1))*deltaTcbe/2;
            sumVy=sumVy+(imuData.accY(j)+imuData.accY(j+1))*deltaTcbe/2;
            sumVz=sumVz+(imuData.accZ(j)+imuData.accZ(j+1))*deltaTcbe/2;
        end
    end
    B4Cbe(iter4Cbe*3+1,1)=sumVx;
    B4Cbe(iter4Cbe*3+1,2)=sumVy;
    B4Cbe(iter4Cbe*3+1,3)=sumVz;
    B4Cbe(iter4Cbe*3+2,4)=sumVx;
    B4Cbe(iter4Cbe*3+2,5)=sumVy;
    B4Cbe(iter4Cbe*3+2,6)=sumVz;
    B4Cbe(iter4Cbe*3+3,7)=sumVx;
    B4Cbe(iter4Cbe*3+3,8)=sumVy;
    B4Cbe(iter4Cbe*3+3,9)=sumVz;
    iter4Cbe=iter4Cbe+1;
end
Ccbe=B4Cbe'*B4Cbe\B4Cbe'*L4Cbe;
Ccbe=reshape(Ccbe,3,3);
cbe=(Ccbe+Ccbe')/2;
% construct equation
%近似值
v0=[nodeSol.Vx(nodeSolSize),nodeSol.Vy(nodeSolSize),nodeSol.Vz(nodeSolSize)]';
ps=[nodeSol.X(nodeSolSize),nodeSol.Y(nodeSolSize),nodeSol.Z(nodeSolSize)]';

prT=imuData.gpsSecond(insEndIndex);
prA=[imuData.accX(insEndIndex),imuData.accY(insEndIndex),imuData.accZ(insEndIndex)-9.8]';
sumJ=zeros(3,1);
sumV=zeros(3,1);
sumM=zeros(3,3);
sumT=zeros(3,1);
%观测方程矩阵
B=zeros((insEndIndex-insStartIndex),15);
L=zeros((insEndIndex-insStartIndex),1);
P=zeros(insEndIndex-insStartIndex,insEndIndex-insStartIndex);
X0=[nodeSol.Vx(nodeSolSize),nodeSol.Vy(nodeSolSize),nodeSol.Vz(nodeSolSize)...
    cbe(1,1),cbe(1,2),cbe(1,3),cbe(2,1),cbe(2,2),cbe(2,3),cbe(3,1),cbe(3,2),cbe(3,3),0,0,0]';
prP=ps;
preV=v0;
iterIndex=1;
sumDeltaT=0;
approxP=[];
approxV=[];
Pi1s=[];
Pi2s=[];
Pi3s=[];
Bp1s=[];
Bp2s=[];
Bp3s=[];
for i=(insEndIndex-1):-1:insStartIndex
    %time interval
    deltaT=imuData.gpsSecond(i)-prT;
    %total time eclpise
    sumDeltaT=sumDeltaT+deltaT;
    %curr accelerometer
    curA=[imuData.accX(i),imuData.accY(i),imuData.accZ(i)-9.8]';
    %curr velocity
    curV=preV+cbe*(curA+prA)/2*deltaT;
    %curr position
    curP=prP+(curV+preV)/2*deltaT;
    
    approxP=[approxP;curP(1),curP(2),curP(3)];
    approxV=[approxV;curV(1),curV(2),curV(3)];
    
    
    sumJ=sumJ+cbe/4*(curA+prA)*deltaT*deltaT+sumV*deltaT;
    sumV=sumV+cbe/2*(curA+prA)*deltaT;
    sumM=sumM+cbe/2*deltaT*deltaT+sumT*deltaT;
    sumT=sumT+cbe*deltaT;
    
    %
    diffP=curP-avgXYZ';
    Pi1=d(3)*diffP(2)-d(2)*diffP(3);
    Pi2=-d(3)*diffP(1)+d(1)*diffP(3);
    Pi3=d(2)*diffP(1)-d(1)*diffP(2);
    Pi1s=[Pi1s,Pi1];
    Pi2s=[Pi2s,Pi2];
    Pi3s=[Pi3s,Pi3];
    
    Bp1=2*(d(2)*Pi3-d(3)*Pi2);
    Bp2=2*(d(3)*Pi1-d(1)*Pi3);
    Bp3=2*(-d(2)*Pi1+d(1)*Pi2);
    Bp1s=[Bp1s,Bp1];
    Bp2s=[Bp2s,Bp2];
    Bp3s=[Bp3s,Bp3];
    B(iterIndex,1)=Bp1*sumDeltaT;
    B(iterIndex,2)=Bp2*sumDeltaT;
    B(iterIndex,3)=Bp3*sumDeltaT;
    B(iterIndex,4)=Bp1*sumJ(1);
    B(iterIndex,5)=Bp1*sumJ(2);
    B(iterIndex,6)=Bp1*sumJ(3);
    B(iterIndex,7)=Bp2*sumJ(1);
    B(iterIndex,8)=Bp2*sumJ(2);
    B(iterIndex,9)=Bp2*sumJ(3);
    B(iterIndex,10)=Bp3*sumJ(1);
    B(iterIndex,11)=Bp3*sumJ(2);
    B(iterIndex,12)=Bp3*sumJ(3);
    B(iterIndex,13)=Bp1*sumM(1,1)+Bp2*sumM(2,1)+Bp3*sumM(3,1);
    B(iterIndex,14)=Bp1*sumM(1,2)+Bp2*sumM(2,2)+Bp3*sumM(3,2);
    B(iterIndex,15)=Bp1*sumM(1,3)+Bp2*sumM(2,3)+Bp3*sumM(3,3);
    L(iterIndex)=Pi1*Pi1+Pi2*Pi2+Pi3*Pi3;
    
%     B(iterIndex*3+1,1)=deltaT;
%     B(iterIndex*3+1,4)=sumJ(1);
%     B(iterIndex*3+1,5)=sumJ(2);
%     B(iterIndex*3+1,6)=sumJ(3);
%     B(iterIndex*3+1,13)=sumM(1,1);
%     B(iterIndex*3+1,14)=sumM(1,2);
%     B(iterIndex*3+1,15)=sumM(1,3);
%     B(iterIndex*3+1,:)=*B(iterIndex*3+1,:);
%     
%     B(iterIndex*3+2,2)=deltaT;
%     B(iterIndex*3+2,7)=sumJ(1);
%     B(iterIndex*3+2,8)=sumJ(2);
%     B(iterIndex*3+2,9)=sumJ(3);
%     B(iterIndex*3+2,13)=sumM(2,1);
%     B(iterIndex*3+2,14)=sumM(2,2);
%     B(iterIndex*3+2,15)=sumM(2,3);
%     B(iterIndex*3+2,:)=(d(3)*Pi1-d(1)*Pi3)*B(iterIndex*3+2,:);
% 
%     
% 
%     B(iterIndex*3+3,3)=deltaT;
%     B(iterIndex*3+3,10)=sumJ(1);
%     B(iterIndex*3+3,11)=sumJ(2);
%     B(iterIndex*3+3,12)=sumJ(3);
%     B(iterIndex*3+3,13)=sumM(3,1);
%     B(iterIndex*3+3,14)=sumM(3,2);
%     B(iterIndex*3+3,15)=sumM(3,3);
%     B(iterIndex*3+3,:)=(-d(2)*Pi1+d(1)*Pi2)*B(iterIndex*3+3,:);
    
     P(iterIndex,iterIndex)=exp(-iterIndex);

    prT=imuData.gpsSecond(i);
    prP=curP;
    prA=curA;
    preV=curV;
    iterIndex=iterIndex+1;
end
figure
plot3(approxP(:,1),approxP(:,2),approxP(:,3))
figure
plot(L)
 L=(B*X0-L);                          
X=(B'*P*B)\B'*P*L;
figure
plot(B)
end
end

