close all
clear

%% 构造GPS数据
load gps2018-9-15.mat;
gpsData.gpsWeek=gps{:,1};
gpsData.gpsSecond=gps{:,2};
gpsData.latitude=gps{:,3};
gpsData.longitude=gps{:,4};
gpsData.height=gps{:,5};

gpsData.sdn=gps{:,'sdn'};
gpsData.sde=gps{:,'sde'};
gpsData.sdu=gps{:,'sdu'};


plot(gpsData.longitude(10:60),gpsData.latitude(10:60),'.');
title('The field trajectory');
xlabel('longitude(deg)');
ylabel('latitude(deg)');
figure
plot(gpsData.gpsSecond,gpsData.sdn,'r',gpsData.gpsSecond,gpsData.sde,'g',gpsData.gpsSecond,gpsData.sdu,'b');
title('the std of gps solution');
xlabel('gps second(s)');
ylabel('std(m)');
legend('North', 'East', 'Height');
%% 构造INS数据
load imu2018-9-15.mat;
[imuNum,~]=size(imu);
toRmv=[];
for i=1:imuNum
    if(any(ismissing(imu(i,[3,5,6,7,8,9,10]))))
        toRmv=[toRmv,i];
    end
end
imu(toRmv,:)=[];
[imuNum,~]=size(imu);
imuData.temp=imu{:,'EKFTemp'};
imuData.gX=imu{:,'EKFRawXGyro'};
imuData.gY=imu{:,'EKFRawYGyro'};
imuData.gZ=imu{:,'EKFRawZGyro'};
imuData.accX=imu{:,'EKFRawXAccel'};
imuData.accY=imu{:,'EKFRawYAccel'};
imuData.accZ=imu{:,'EKFRawZAccel'};
imuData.gpsWeek=zeros(imuNum,1);
imuData.gpsSecond=zeros(imuNum,1);
for i=1:imuNum
    if(~ismissing(imu(i,'GPStime')))
        ts=strsplit(string(imu{i,'GPStime'}),':');
        imuData.gpsWeek(i)=str2double(convertStringsToChars(ts(1)));
        imuData.gpsSecond(i)=str2double(convertStringsToChars(ts(2)));
    else 
        imuData.gpsWeek(i)=NaN;
        imuData.gpsSecond(i)=NaN;
    end   
end
figure
plot(imuData.temp)
title('temperature')
figure 
subplot(3,1,1)
plot(imuData.accX)
title('accX')
subplot(3,1,2)
plot(imuData.accY)
title('accY')
subplot(3,1,3)
plot(imuData.accZ)
title('accZ')

figure 
subplot(3,1,1)
plot(imuData.gX)
title('gX')
subplot(3,1,2)
plot(imuData.gY)
title('gY')
subplot(3,1,3)
plot(imuData.gZ)
title('gZ')
%%