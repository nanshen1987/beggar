function [pass,bias,distance,cbe] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins)
%% check condition for allignment
%   nodeSol gpsNode for check
%   imuData ins data
%   gpsWinsize allgin window
%   insGpsRatio ratio  of gps and ins
%   epoch4ins current epoch of ins

limitBias=10;
limitDistance=15;
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
    cbe=eye(3);
    return;
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
cbeM=zeros(3,3);
cbeM(1,:)=Ccbe(1:3);
cbeM(2,:)=Ccbe(4:6);
cbeM(3,:)=Ccbe(7:9);
cbe=cbeM;


%观测方程矩阵
%2018.10.23变更，只求姿态，忽略初始位置和初始速度（不作为待估值求解）
B=zeros((insEndIndex-insStartIndex),9);
L=zeros((insEndIndex-insStartIndex),1);
iterIndex=1;
for i=(insEndIndex-1):-1:insStartIndex
    %找出离该点最近的GPS点
    insSec=imuData.gpsSecond(i);
    deltaInsNode=1000000;
    matchNodeIndex=0;
    matchNodeSec=0;
    for j=nodeSolSize-gpsWinsize+1:nodeSolSize
        nodeSec=nodeSol.gpsSecond(j);
        if abs(nodeSec-insSec)<deltaInsNode
            deltaInsNode=abs(nodeSec-insSec);
            matchNodeIndex=j;
            matchNodeSec=nodeSec;
        end
    end
    %估计该点坐标
    %判断GPS点靠前还是靠后
    if matchNodeSec>=insSec
        i_end=i;
        for k=i:insEndIndex
            if imuData.gpsSecond(k)<=matchNodeSec
                i_end=k;
            end
        end
        p0x=nodeSol.X(matchNodeIndex);
        p0y=nodeSol.Y(matchNodeIndex);
        p0z=nodeSol.Z(matchNodeIndex);
        v0x=nodeSol.Vx(matchNodeIndex);
        v0y=nodeSol.Vy(matchNodeIndex);
        v0z=nodeSol.Vz(matchNodeIndex);
        preT=matchNodeSec;
        preP=[p0x p0y p0z]';
        preV=[v0x v0y v0z]';
        preA=[imuData.accX(i_end) imuData.accY(i_end) imuData.accZ(i_end)]';
        sumJ=zeros(3,1);
        sumV=zeros(3,1);
        for k=i_end:i
            deltaT=imuData.gpsSecond(k)-preT;
            curA=[imuData.accX(k) imuData.accY(k) imuData.accZ(k)]';
            curV=preV+cbe*(preA+curA)/2*deltaT;
            curP=preP+(preV+curV)/2.0*deltaT;
            sumJ=sumJ+cbe/4*(curA+preA)*deltaT*deltaT+sumV*deltaT;
            sumV=sumV+cbe/2*(curA+preA)*deltaT;
            
            preT=imuData.gpsSecond(k);
            preP=curP;
            preV=curV;
            preA=curA;
        end
        
        diffP=preP-avgXYZ';
        Pi1=d(3)*diffP(2)-d(2)*diffP(3);
        Pi2=-d(3)*diffP(1)+d(1)*diffP(3);
        Pi3=d(2)*diffP(1)-d(1)*diffP(2);
        Bp1=2*(d(2)*Pi3-d(3)*Pi2);
        Bp2=2*(d(3)*Pi1-d(1)*Pi3);
        Bp3=2*(-d(2)*Pi1+d(1)*Pi2);
        
        B(iterIndex,1)=Bp1*sumJ(1);
        B(iterIndex,2)=Bp1*sumJ(2);
        B(iterIndex,3)=Bp1*sumJ(3);
        B(iterIndex,4)=Bp2*sumJ(1);
        B(iterIndex,5)=Bp2*sumJ(2);
        B(iterIndex,6)=Bp2*sumJ(3);
        B(iterIndex,7)=Bp3*sumJ(1);
        B(iterIndex,8)=Bp3*sumJ(2);
        B(iterIndex,9)=Bp3*sumJ(3);
        L(iterIndex)=Pi1*Pi1+Pi2*Pi2+Pi3*Pi3;
    else
         i_start=i;
        for k=i:-1:insStartIndex
            if imuData.gpsSecond(k)>matchNodeSec
                i_start=k;
            end
        end
        p0x=nodeSol.X(matchNodeIndex);
        p0y=nodeSol.Y(matchNodeIndex);
        p0z=nodeSol.Z(matchNodeIndex);
        v0x=nodeSol.Vx(matchNodeIndex);
        v0y=nodeSol.Vy(matchNodeIndex);
        v0z=nodeSol.Vz(matchNodeIndex);
        preT=matchNodeSec;
        preP=[p0x p0y p0z]';
        preV=[v0x v0y v0z]';
        preA=[imuData.accX(i_start) imuData.accY(i_start) imuData.accZ(i_start)]';
        sumJ=zeros(3,1);
        sumV=zeros(3,1);
        for k=i_start:i
            deltaT=imuData.gpsSecond(k)-preT;
            curA=[imuData.accX(k) imuData.accY(k) imuData.accZ(k)]';
            curV=preV+cbe*(preA+curA)/2*deltaT;
            curP=preP+(preV+curV)/2.0*deltaT;
            sumJ=sumJ+cbe/4*(curA+preA)*deltaT*deltaT+sumV*deltaT;
            sumV=sumV+cbe/2*(curA+preA)*deltaT;
            
            preT=imuData.gpsSecond(k);
            preP=curP;
            preV=curV;
            preA=curA;
        end
        
        diffP=preP-avgXYZ';
        Pi1=d(3)*diffP(2)-d(2)*diffP(3);
        Pi2=-d(3)*diffP(1)+d(1)*diffP(3);
        Pi3=d(2)*diffP(1)-d(1)*diffP(2);
        Bp1=2*(d(2)*Pi3-d(3)*Pi2);
        Bp2=2*(d(3)*Pi1-d(1)*Pi3);
        Bp3=2*(-d(2)*Pi1+d(1)*Pi2);
        
        B(iterIndex,1)=Bp1*sumJ(1);
        B(iterIndex,2)=Bp1*sumJ(2);
        B(iterIndex,3)=Bp1*sumJ(3);
        B(iterIndex,4)=Bp2*sumJ(1);
        B(iterIndex,5)=Bp2*sumJ(2);
        B(iterIndex,6)=Bp2*sumJ(3);
        B(iterIndex,7)=Bp3*sumJ(1);
        B(iterIndex,8)=Bp3*sumJ(2);
        B(iterIndex,9)=Bp3*sumJ(3);
        L(iterIndex)=Pi1*Pi1+Pi2*Pi2+Pi3*Pi3;
    end
    iterIndex=iterIndex+1;
end
% figure
% plot(B)
% X=(B'*B)\B'*L;
% cbe=cbe+reshape(X,3,3)';
end
end

