close all
clear
%% ����GPS����
load gps2018-10-8.mat;
gpsData.gpsWeek=gps{:,1};
gpsData.gpsSecond=gps{:,2};
gpsData.X=gps{:,3};
gpsData.Y=gps{:,4};
gpsData.Z=gps{:,5};

gpsData.sdx=gps{:,'sdx'};
gpsData.sdy=gps{:,'sdy'};
gpsData.sdz=gps{:,'sdz'};

gpsData.sdxy=gps{:,'sdxy'};
gpsData.sdyz=gps{:,'sdyz'};
gpsData.sdzx=gps{:,'sdzx'};

gpsData.Vx=gps{:,'Vx'};
gpsData.Vy=gps{:,'Vy'};
gpsData.Vz=gps{:,'Vz'};


scatter3(gpsData.X,gpsData.Y,gpsData.Z,'.');
title('The field trajectory');
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');

figure
plot(gpsData.gpsSecond,gpsData.Vx,'r',gpsData.gpsSecond,gpsData.Vy,'g',gpsData.gpsSecond,gpsData.Vz,'b');
title('the velocity of gps solution');
xlabel('gps second(s)');
ylabel('velocity(m/s)');
legend('X', 'Y', 'Z');
[gpsNum,~]=size(gpsData.gpsSecond);
%% ����INS����
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
subplot(3,1,2)
plot(imuData.gpsSecond,imuData.accY);
subplot(3,1,3)
plot(imuData.gpsSecond,imuData.accZ);
figure 
subplot(3,1,1)
plot(imuData.gpsSecond,imuData.gX);
subplot(3,1,2)
plot(imuData.gpsSecond,imuData.gY);
subplot(3,1,3)
plot(imuData.gpsSecond,imuData.gZ);


%% self allignment and filter

gpsWinsize=10;
insGpsRatio=2;
minDistance=50;
minInterval=3;
%���н�
sol.X=[];
sol.Y=[];
sol.Z=[];
sol.Vx=[];
sol.Vy=[];
sol.Vz=[];
sol.type=[];%0 gps 1 INS 2 INS+GPS
sol.gpsSecond=[];
%���߽ڵ�� type 0 gps  2 INS+GPS
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
fixed=0;%��¼ϵͳ��ǰ�̶�״̬
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
            [pass,bias,distance] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins);
            biases=[biases,bias];
            distances=[distances,distance];      
            fixed=pass          
        end
    else
        nodeSol.X=[nodeSol.X,gpsData.X(epoch4gps)];
        nodeSol.Y=[nodeSol.Y,gpsData.Y(epoch4gps)];
        nodeSol.Z=[nodeSol.Z,gpsData.Z(epoch4gps)];
        nodeSol.Vx=[nodeSol.Vx,gpsData.Vx(epoch4gps)];
        nodeSol.Vy=[nodeSol.Vy,gpsData.Vy(epoch4gps)];
        nodeSol.Vz=[nodeSol.Vz,gpsData.Vz(epoch4gps)];
        nodeSol.type=[nodeSol.type,0];
        nodeSol.gpsSecond=[nodeSol.gpsSecond,gpsData.gpsSecond(epoch4gps)];
        [pass,bias,distance] = check4Allign(nodeSol,imuData,gpsWinsize,insGpsRatio,epoch4ins);
        %%
        if pass==1 && fixed==1
        %ins
            
        elseif pass==1 && fixed==0
        %gps
                
            fixed=1;
        elseif pass==0 && fixed==1
         %gps ins   
            fixed=0;        
        end
        biases=[biases,bias];
        distances=[distances,distance];
        
    end
end
figure
mean(biases)
subplot(3,1,1);
plot(biases,'r');
title('fit bias')
mean(distances)
subplot(3,1,2);
plot(distances,'g');
title('chnage distance')
subplot(3,1,3);
mean(distances./biases)
plot(distances./biases,'b')
title('distances./biases')


