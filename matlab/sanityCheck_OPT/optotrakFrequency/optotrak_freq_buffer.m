% I changed the mdule to 3020. 
% Use the buffer mode to collect data, pin 4 is where the buffer was
% recorded. 
% The sampling rate was using 200Hz and the sampled time recorded in
% TS.mat, whereas the data recorded in the intermediate data (by buffered) 

% cpDatarg2(4121)

% See the pulse time here 
clear; 
ss_num = 4124;
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
    'KingKongTSync.0' num2str(ss_num) '.mat']);

time_idx = data.eventsL == 5; % OPT
trial_idx = data.eventsTrials(time_idx);
time = data.eventsT(time_idx);

% compare the optotrak buffer data  
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/' ...
    'KingKong.0' num2str(ss_num) '.mat']);
databuffer = Data.QL.Data.OPTO_BUFFER_DATA;

trials_list = unique(trial_idx);

for triali = 1:3
% for triali = 2:8 
    trial = trials_list(triali);
    num_datapoints = (length(databuffer{triali})/8 - 1)/3  % # of buffers
    num_pulses = sum(trial_idx == trial);            % # of frames 
    ratio = num_datapoints - num_pulses;            % # diff
end

%% convert the data to float 
data_uint8 = Data.QL.Data.OPTO_BUFFER_DATA{2};


dat_arr = zeros(1,length(data_uint8)/8);
for byte_i = 1:length(data_uint8)/8
    uint8_arr = [];
    offset = (byte_i-1)*8;
    for dat_pos_i = 1:8
    %    uint8_arr = [uint8_arr dec2bin(typecast(int8(data_uint8(dat_pos_i)),'uint8'),8)]
        uint8_arr = [dec2bin(typecast(int8(data_uint8(offset+dat_pos_i)),'uint8'),8) uint8_arr];
    end
    q = quantizer('double');
    B = bin2num(q, uint8_arr);
    dat_arr(byte_i) = B;
end
pos = reshape(dat_arr(2:end), [3, length(dat_arr(2:end))/3]);
plot(pos')'

for dat_pos_i = 1:8
    typecast(int8(data_uint8(dat_pos_i)),'uint8')

end


