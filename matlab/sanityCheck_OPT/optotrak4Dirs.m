% 2022-10-03 
% after changing optotrak mounting and the subject-to-robot configuration, 
% I'm testing if optotrak works 

cpDatarg2(4377)

sstmp = SessionScan(4377); 

figure(); 
axh(1) = subplot(4,1,1);
plot(sstmp.data.t, sstmp.data.f(1,:));
axh(2) = subplot(4,1,2);
% plot(sstmp.data.t, sstmp.data.optx(1:3,:) - mean(sstmp.data.optx(1:3,:), 'omitnan'));
plot(sstmp.opt.datah.t, sstmp.opt.datah.x(1:3,:) - mean(sstmp.opt.datah.x(1:3,:),2, 'omitnan'));
axh(3) = subplot(4,1,3);
axh(4) = subplot(4,1,4);