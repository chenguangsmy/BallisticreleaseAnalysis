% Himanshu found wired noise in his ATI FT, I'm checking if I has the same
% thing

% ss4357, raw data
% ss4358, raw data with online Lowpass filter 3000Hz

ftmp1 = SessionScanFT(4357);
ftmp2 = SessionScanFT(4358);
%%
figure('name', 'force'); 
for xyz_i = 1:3
    axh(xyz_i,1) = subplot(3,2,(xyz_i-1)*2+1); 
%     plot(ftmp1.elapse - ftmp1.elapse(1), ftmp1.force_origin(xyz_i,:));
    plot(ftmp1.elapse - ftmp1.elapse(1), ftmp1.force_origin(xyz_i,:) - ftmp1.force_origin(xyz_i,end));
    title('raw force');
    xlabel('time from start (s)');
    ylabel('Force (N)');
end
for xyz_i = 1:3
    axh(xyz_i,2) = subplot(3,2,(xyz_i-1)*2+2);  
%     plot(ftmp2.elapse - ftmp2.elapse(1), ftmp2.force_origin);
    plot(ftmp2.elapse - ftmp2.elapse(1), ftmp2.force_origin(xyz_i,:) - ftmp2.force_origin(xyz_i,end));
    title('online filtered force 3000Hz');
    xlabel('time from start (s)');
    ylabel('Force (N)');
end
linkaxes(axh(:), 'xy');

figure('name', 'torque'); 
for xyz_i = 1:3
    axh(xyz_i,1) = subplot(3,2,(xyz_i-1)*2+1); 
%     plot(ftmp1.elapse - ftmp1.elapse(1), ftmp1.torque_origin(xyz_i,:));
    plot(ftmp1.elapse - ftmp1.elapse(1), ftmp1.torque_origin(xyz_i,:) - ftmp1.torque_origin(xyz_i,end));
    title('raw force');
    xlabel('time from start (s)');
    ylabel('T (N/m)');
end
for xyz_i = 1:3
    axh(xyz_i,2) = subplot(3,2,(xyz_i-1)*2+2);  
%     plot(ftmp2.elapse - ftmp2.elapse(1), ftmp2.torque_origin(xyz_i,:));
    plot(ftmp2.elapse - ftmp2.elapse(1), ftmp2.torque_origin(xyz_i,:) - ftmp2.torque_origin(xyz_i,end));
    title('online filtered force 3000Hz');
    xlabel('time from start (s)');
    ylabel('T (N/m)');
end
linkaxes(axh(:), 'xy');