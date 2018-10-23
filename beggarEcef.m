close all
clear
clc
%% 构造GPS数据
load gps2018-10-8.mat;
smapled=1:1:800;
gpsData.gpsWeek=gps{smapled,1};
gpsData.gpsSecond=gps{smapled,2};
gpsData.X=gps{smapled,3};
gpsData.Y=gps{smapled,4};
gpsData.Z=gps{smapled,5};

gpsData.sdx=gps{smapled,'sdx'};
gpsData.sdy=gps{smapled,'sdy'};
gpsData.sdz=gps{smapled,'sdz'};

gpsData.sdxy=gps{smapled,'sdxy'};
gpsData.sdyz=gps{smapled,'sdyz'};
gpsData.sdzx=gps{smapled,'sdzx'};

gpsData.Vx=gps{smapled,'Vx'};
gpsData.Vy=gps{smapled,'Vy'};
gpsData.Vz=gps{smapled,'Vz'};



figure
plot(gpsData.gpsSecond,gpsData.Vx,'r',gpsData.gpsSecond,gpsData.Vy,'g',gpsData.gpsSecond,gpsData.Vz,'b');
title('the velocity of gps solution');
xlabel('gps second(s)');
ylabel('velocity(m/s)');
legend('X', 'Y', 'Z');
[gpsNum,~]=size(gpsData.gpsSecond);
%% 构造INS数据
load imu2018-10-8.mat;
[imuRawNum,~]=size(imu);

imuData.gX=zeros(imuRawNum,1);
imuData.gY=zeros(imuRawNum,1);
imuData.gZ=zeros(imuRawNum,1);
imuData.accX=zeros(imuRawNum,1);
imuData.accY=zeros(imuRawNum,1);
imuData.accZ=zeros(imuRawNum,1);
imuData.iTOW=zeros(imuRawNum,1);
imuData.gpsWeek=zeros(imuRawNum,1);
imuData.gpsSecond=zeros(imuRawNum,1);
imuNum=0;
for i=1:imuRawNum
    if(any(ismissing(imu(i,[3,4,5,6,7,8]))))
       continue;
    end
    imuNum=imuNum+1;
    imuData.gX(imuNum)=imu{i,'EKFRawXGyro'};
    imuData.gY(imuNum)=imu{i,'EKFRawYGyro'};
    imuData.gZ(imuNum)=imu{i,'EKFRawZGyro'};
    imuData.accX(imuNum)=imu{i,'EKFRawXAccel'};
    imuData.accY(imuNum)=imu{i,'EKFRawYAccel'};
    imuData.accZ(imuNum)=imu{i,'EKFRawZAccel'};
    imuData.iTOW(imuNum)=imu{i,'iTOW'};
    if(string(imu{i,'GPStime'})=="")
        %foreward  fix it
        imuData.gpsWeek(imuNum)=imuData.gpsWeek(imuNum-1);
        imuData.gpsSecond(imuNum)=imuData.gpsSecond(imuNum-1)+imuData.iTOW(imuNum)-imuData.iTOW(imuNum-1);
    else 
        ts=strsplit(string(imu{i,'GPStime'}),':');
        imuData.gpsWeek(imuNum)=str2double(convertStringsToChars(ts(1)));
        imuData.gpsSecond(imuNum)=str2double(convertStringsToChars(ts(2)));
    end   
end
imuData.gX=imuData.gX(1:imuNum);
imuData.gY=imuData.gY(1:imuNum);
imuData.gZ=imuData.gZ(1:imuNum);
imuData.accX=imuData.accX(1:imuNum);
imuData.accY=imuData.accY(1:imuNum);
imuData.accZ=imuData.accZ(1:imuNum);
imuData.iTOW=imuData.iTOW(1:imuNum);
imuData.gpsWeek=imuData.gpsWeek(1:imuNum);
imuData.gpsSecond=imuData.gpsSecond(1:imuNum);
figure 
subplot(3,1,1)
plot(imuData.gpsSecond,imuData.accX);
title('imuData.accX')
subplot(3,1,2)
plot(imuData.gpsSecond,imuData.accY);
title('imuData.accY')
subplot(3,1,3)
plot(imuData.gpsSecond,imuData.accZ);
title('imuData.accZ')
figure 
subplot(3,1,1)
plot(imuData.gpsSecond,imuData.gX);
title('imuData.gX')
subplot(3,1,2)
plot(imuData.gpsSecond,imuData.gY);
title('imuData.gY')
subplot(3,1,3)
plot(imuData.gpsSecond,imuData.gZ);
title('imuData.gZ')


%% self allignment and filter

gpsWinsize=8;
insGpsRatio=2;
minDistance=5;
minInterval=3;
%所有解
sol.X=[];
sol.Y=[];
sol.Z=[];
sol.Vx=[];
sol.Vy=[];
sol.Vz=[];
sol.type=[];%0 gps 1 INS 2 INS+GPS
sol.gpsSecond=[];
%决策节点解 type 0 gps  2 INS+GPS
nodeSol.X=[];
nodeSol.Y=[];
nodeSol.Z=[];
nodeSol.Vx=[];
nodeSol.Vy=[];
nodeSol.Vz=[];
nodeSol.type=[];%0 gps  2 INS+GPS
nodeSol.gpsSecond=[];
epoch4ins=1;
biases=[];
distances=[];
fixed=0;%记录系统当前固定状态
cbe=eye(3);
for epoch4gps=1:gpsNum
    %
    if epoch4gps<=gpsWinsize
        sol.X=[sol.X,gpsData.X(epoch4gps)];
        sol.Y=[sol.Y,gpsData.Y(epoch4gps)];
        sol.Z=[sol.Z,gpsData.Z(epoch4gps)];
        sol.Vx=[sol.Vx,gpsData.Vx(epoch4gps)];
        sol.Vy=[sol.Vy,gpsData.Vy(epoch4gps)];
        sol.Vz=[sol.Vz,gpsData.Vz(epoch4gps)];
        sol.type=[sol.type,0];
        sol.gpsSecond=[sol.gpsSecond,gpsData.gpsSecond(epoch4gps)];
        
        nodeSol.X=[nodeSol.X,gpsData.X(epoch4gps)];
        nodeSol.Y=[nodeSol.Y,gpsData.Y(epoch4gps)];
        nodeSol.Z=[nodeSol.Z,gpsData.Z(epoch4gps)];
        nodeSol.Vx=[nodeSol.Vx,gpsData.Vx(epoch4gps)];
        nodeSol.Vy=[nodeSol.Vy,gpsData.Vy(epoch4gps)];
        nodeSol.Vz=[nodeSol.Vz,gpsData.Vz(epoch4gps)];
        nodeSol.type=[nodeSol.type,0];
        nodeSol.gpsSecond=[nodeSol.gpsSecond,gpsData.gpsSecond(epoch4gps)];
        if epoch4gps==gpsWinsize
            [pass,bias,distance,tcbe] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins);
            biases=[biases,bias];
            distances=[distances,distance];      
            fixed=pass 
            if pass
                cbe=tcbe;
            end
        end
    else
        sol.X=[sol.X,gpsData.X(epoch4gps)];
        sol.Y=[sol.Y,gpsData.Y(epoch4gps)];
        sol.Z=[sol.Z,gpsData.Z(epoch4gps)];
        sol.Vx=[sol.Vx,gpsData.Vx(epoch4gps)];
        sol.Vy=[sol.Vy,gpsData.Vy(epoch4gps)];
        sol.Vz=[sol.Vz,gpsData.Vz(epoch4gps)];
        sol.type=[sol.type,0];
        sol.gpsSecond=[sol.gpsSecond,gpsData.gpsSecond(epoch4gps)];
        
        
        nodeSol.X=[nodeSol.X,gpsData.X(epoch4gps)];
        nodeSol.Y=[nodeSol.Y,gpsData.Y(epoch4gps)];
        nodeSol.Z=[nodeSol.Z,gpsData.Z(epoch4gps)];
        nodeSol.Vx=[nodeSol.Vx,gpsData.Vx(epoch4gps)];
        nodeSol.Vy=[nodeSol.Vy,gpsData.Vy(epoch4gps)];
        nodeSol.Vz=[nodeSol.Vz,gpsData.Vz(epoch4gps)];
        nodeSol.type=[nodeSol.type,0];
        nodeSol.gpsSecond=[nodeSol.gpsSecond,gpsData.gpsSecond(epoch4gps)];
        [pass,bias,distance,tcbe] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins);
       
        if pass==1 && fixed==1
        %ins
            cbe=tcbe;
        elseif pass==1 && fixed==0
        %gps
            cbe=tcbe; 
            fixed=1;
        elseif pass==0 && fixed==1
         %gps ins   
            fixed=0;        
        end
        biases=[biases,bias];
        distances=[distances,distance];
        
    end
    %% check if cbe is fixed for ins
    if fixed&&epoch4gps<gpsNum
        curGpsSec=gpsData.gpsSecond(epoch4gps);
        nextGpsSec=gpsData.gpsSecond(epoch4gps+1);
        preP=[gpsData.X(epoch4gps) gpsData.Y(epoch4gps) gpsData.Z(epoch4gps)]';
        preV=[gpsData.Vx(epoch4gps) gpsData.Vy(epoch4gps) gpsData.Vz(epoch4gps)]';
        preT=curGpsSec;
        preA=[];
        for j=epoch4ins:imuNum
            curInsSec=imuData.gpsSecond(j);
            if curInsSec>curGpsSec && curInsSec<nextGpsSec
                deltaT=curInsSec-preT;
                curA=[imuData.accX(j),imuData.accY(j),imuData.accZ(j)]';
                if isempty(preA)
                    curV=preV+cbe*curA*deltaT;
                else
                    curV=preV+cbe*(curA+preA)/2.0*deltaT;
                end
                curP=preP+(preV+curV)/2.0*deltaT;
                                
                sol.X=[sol.X,curP(1)];
                sol.Y=[sol.Y,curP(2)];
                sol.Z=[sol.Z,curP(3)];
                sol.Vx=[sol.Vx,curV(1)];
                sol.Vy=[sol.Vy,curV(2)];
                sol.Vz=[sol.Vz,curV(3)];
                sol.type=[sol.type,1];
                sol.gpsSecond=[sol.gpsSecond,curInsSec];
                
                preA=curA;
                preV=curV;
                preP=curP;
                preT=curInsSec;
            end
        end
    end
end
% figure
% mean(biases)
% subplot(3,1,1);
% plot(biases,'r');
% title('fit bias')
% mean(distances)
% subplot(3,1,2);
% plot(distances,'g');
% title('chnage distance')
% subplot(3,1,3);
% mean(distances./biases)
% plot(distances./biases,'b')
% title('distances./biases')
figure
scatter3(gpsData.X,gpsData.Y,gpsData.Z,'g.');
title('The field trajectory');
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
hold on 
scatter3(sol.X,sol.Y,sol.Z,'b+');
figure
plot(gpsData.gpsSecond,gpsData.X,'go');
hold on 
plot(sol.gpsSecond,sol.X,'b+');
title('PK X');
figure
plot(gpsData.gpsSecond,gpsData.Y,'go');
hold on 
plot(sol.gpsSecond,sol.Y,'b+');
title('PK Y');
figure
plot(gpsData.gpsSecond,gpsData.Z,'go');
hold on 
plot(sol.gpsSecond,sol.Z,'b+');
title('PK Z');


