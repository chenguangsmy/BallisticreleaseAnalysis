%% Test non-linearity of robot using spring test 
% 1. No-calibrated. 
ss2333 = SessionScan(2333); 
ss2336 = SessionScan(2336);
axhv1 = ss2333.plotMeantrialVel_sameCond_overlap(1, 1, 1);
axhv1 = ss2336.plotMeantrialVel_sameCond_overlap(-1, axhv1, 2);
ylim([-0.7 0.7]);
title('movement difference before calibration');
legend('forward', 'backward (invert)');
% 2. After-calibration
ss2338 = SessionScan(2338);
ss2340 = SessionScan(2340);
axhv2 = ss2340.plotMeantrialVel_sameCond_overlap(1, 1, 2);
axhv2 = ss2338.plotMeantrialVel_sameCond_overlap(-1, axhv2, 1);
ylim([-0.7 0.7]);
title('movement difference after calibration');
legend('forward', 'backward (invert)');
%% 3. New sessions
ssTestList = [2341, 2342];
ssFront = SessionScan(ssTestList(1));
ssBack = SessionScan(ssTestList(2));
axhv2 = ssFront.plotMeantrialVel_sameCond_overlap(1, 1, 1);
axhv2 = ssBack.plotMeantrialVel_sameCond_overlap(-1, axhv2, 2);
ylim([-0.7 0.7]);
title('movement difference after calibration');
legend('forward', 'backward (invert)');
