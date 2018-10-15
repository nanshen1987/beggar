function [pass,bias,distance] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins)
%% check condition for allignment
%   nodeSol gpsNode for check
%   imuData ins data
%   gpsWinsize allgin window
%   insGpsRatio ratio  of gps and ins
%   epoch4ins current epoch of ins

limitBias=9;
limitDistance=10;
limitR=2.8641;

[~,nodeSolSize]=size(nodeSol.X);
[~,imuDataSize]=size(imuData.gpsSecond);
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
        if imuData.gpsSecond(i)<=gpsEndSec
            insEndIndex=i;
        else
            break;
        end
    end
    for i=epoch4ins:-1:1
        if imuData.gpsSecond(i)>=gpsStartSec
            insStartIndex=i;
        else
            break;
        end
    end
end
% construct equation
%½üËÆÖµ
v0=[nodeSol.Vx(nodeSolSize),nodeSol.Vy(nodeSolSize),nodeSol.Vz(nodeSolSize)]';
ps=[nodeSol.X(nodeSolSize),nodeSol.Y(nodeSolSize),nodeSol.Z(nodeSolSize)]';
cbe=eyes(3);
prT=imuData.gpsSecond(insEndIndex);
prA=[imuData.accX(insEndIndex),imuData.accY(insEndIndex),imuData.accZ(insEndIndex)]';
sumJ=zeros(3,1);
sumM=zeros(3,3);
prP=ps;
preV=v0;
for i=insEndIndex:-1:insStartIndex
    deltaT=imuData.gpsSecond(insEndIndex)-prT;
    curA=[imuData.accX(i),imuData.accY(i),imuData.accZ(i)]';
    curV=preV-cbe*(curA+prA)/2*deltaT;
    curP=prP-(curV+preV)/2*deltaT;
    
    sumJ=sumJ+(prA+curA)*deltaT;
    
    prT=imuData.gpsSecond(insEndIndex);
    prP=curP;
end
end

