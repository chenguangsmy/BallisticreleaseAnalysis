InterKK = load('../KingKongData/Intermediate/KingKong.01217.mat')

close all
figure
d1 = 2006;
d2=696;
plot(1:d1,InterKK.Data.QL.Data.FORCE_SENSOR_DATA.data(1,:));
hold on;
plot(1:d1,InterKK.Data.QL.Data.FORCE_SENSOR_DATA.data(2,:),'r');
plot(1:d1,InterKK.Data.QL.Data.FORCE_SENSOR_DATA.data(3,:),'g');

figure
plot(1:d2,InterKK.Data.QL.Data.BURT_STATUS.pos_x(1,:));
hold on;
plot(1:d2,InterKK.Data.QL.Data.BURT_STATUS.pos_y(1,:),'r');
plot(1:d2,InterKK.Data.QL.Data.BURT_STATUS.pos_z(1,:),'g');

figure;
plot3(InterKK.Data.QL.Data.BURT_STATUS.pos_x(1, 1:120) + 0.12,...
    InterKK.Data.QL.Data.BURT_STATUS.pos_y(1, 1:120) - 0.5,...
    InterKK.Data.QL.Data.BURT_STATUS.pos_z(1, 1:120) - 0.25);
ylim([-0.03, 0.03]);
xlim([-0.03, 0.03]);
zlim([-0.03,0.03]);
figure; 
plot3(InterKK.Data.QL.Data.BURT_STATUS.pos_x(1, 1:120),...
    InterKK.Data.QL.Data.BURT_STATUS.pos_y(1, 1:120),...
    InterKK.Data.QL.Data.BURT_STATUS.pos_z(1, 1:120));