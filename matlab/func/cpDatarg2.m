function [cp_flag] = cpDatarg2(ss_num)
%COPOYDATARG2 copy data from rg2 to my current computer with fixed source and desiny location 
%   srcDir: '/Volumes/rg2/data/KingKong/...'
%   dstDir: '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data'
%   copy may have one or multiple files. 

srcdir1 = '/Volumes/rg2/data/KingKong/Raw/';
srcdir2 = '/Volumes/rg2/data/KingKong/Formatted/';
srcdir3 = '/Volumes/rg2/data/KingKong/Intermediate/';
srcdirvd = '/Volumes/rg2/data/VideoLoggerData/';

dstdir1 = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/';
dstdir2 = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/';
dstdirvd= '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/VideoLoggerData/';

ss_num = ss_num(:);
cp_flag = zeros(size(ss_num));
% Formatted file, intermediate file, force transducer, WAM, TimeSync, EMG,
% Optotrak, Videologger
% FMT  INT  FT  WAM  TS  EMG  OPT  VDL
cp_flag_str = '11111111';

ERRCODE_LOSTFORMAT = 8;
ERRCODE_LOSTIMTERM = 4;
ERRCODE_LOSTFTCSV  = 2;
ERRCODE_LOSTWAMCSV = 1;
ERRCODE_LOSTTIMSYN = 16;


for ss_i = 1:length(ss_num)
    sprintf('copying... ');
    sprintf('%d / %d  ', ss_i, length(ss_num));
    
    % formatted data
    fname3 = sprintf([srcdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    if (~exist(fname3, 'file'))
%         disp(['no formatted *.mat file for ss' num2str(ss_num(ss_i))]);
        cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTFORMAT;
        cp_flag_str(1) = '0';
    end
    
    % intermediate data
    fname4 = sprintf([srcdir3 'KingKong.%05d.mat'], ss_num(ss_i));
    if (~exist(fname4, 'file'))
%         disp(['no intermediate *.mat file for ss' num2str(ss_num(ss_i))]);
        cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTIMTERM;
        cp_flag_str(2) = '0';
    end
    
    % FT data
    fname1 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname10= sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT0.csv'], ss_num(ss_i));
    if (~exist(fname1, 'file'))
%         disp(['no processed FT file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname10, 'file'))
%             disp('cautiion: no FT file at all');
            cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTFTCSV;
            cp_flag_str(3) = '0';
        else
%             disp('use FT0.csv file.');
        end
    end
    
    % WAM data
    fname2 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname20= sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM0.csv'], ss_num(ss_i));
    fname21= sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAMbin'], ss_num(ss_i));
    if (~exist(fname2, 'file'))
%         disp(['no processed WAM file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname20, 'file'))
%             disp('cautiion: no WAM file at all');
            cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTWAMCSV;
            cp_flag_str(4) = '0';
        else
%             disp('use WAM0.csv file.');
            fname2 = fname20;
        end
        if (exist(fname21, 'file')) % binary file 
            disp('cautiion: WAM file is binary');
        end
    end
    
    % timeSync data
    fname5 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongTimeSync0.mat'], ss_num(ss_i));
    fname51= sprintf([srcdir1 'KingKong.DK.%05d/KingKongTimeSync%05d.mat'], ss_num(ss_i), ss_num(ss_i));
    if (~exist(fname5, 'file'))
%         disp(['no TimeSync *.mat file for ss' num2str(ss_num(ss_i))]);
        fname5 = fname51;
        if (~exist(fname5, 'file'))
            cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTTIMSYN;
            cp_flag_str(5) = '0';
        end
    end
    
    % EMG data
    fname6 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongEMG%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname60= sprintf([srcdir1 'KingKong.DK.%05d/KingKongEMG0.csv'], ss_num(ss_i));
    fname61= sprintf([srcdir1 'KingKong.DK.%05d/KingKongEMG%05d.mat'], ss_num(ss_i), ss_num(ss_i));
    if (~exist(fname6, 'file'))
%         disp(['no processed EMG file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname60, 'file'))
%             disp('cautiion: no EMG file at all');
            if (~exist(fname61, 'file'))
                cp_flag_str(6) = '0';
            end
%             cp_flag(ss_i) = cp_flag(ss_i) + ERRORCODE_LOSTWAMCSV;
        else
%             disp('use EMG0.csv file.');
        end
    end
    
    % OPTOTRAK data
    fname7 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongOPT%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname70= sprintf([srcdir1 'KingKong.DK.%05d/KingKongOPT0.csv'], ss_num(ss_i));
    if (~exist(fname7, 'file'))
%         disp(['no processed OPTOTRAK file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname70, 'file'))
%             disp('cautiion: no OPTOTRAK file at all');
            cp_flag_str(7) = '0';
%             cp_flag(ss_i) = cp_flag(ss_i) + ERRORCODE_LOSTWAMCSV;
        else
            disp('use OPT0.csv file.');
        end
    end
    
    % Videologger data
    fnamevd = sprintf([srcdirvd 'KingKong.DK.' num2str(ss_num(ss_i)) '.mp4']); % video
    fnamevdd = sprintf([srcdirvd 'KingKong.DK.' num2str(ss_num(ss_i)) '.txt']); % video description 
    if (~exist(fnamevd,'file'))
        cp_flag_str(8) = '0';
    end
    
    
    
    dstn1  = sprintf([dstdir1 'KingKongFT%05d.csv'], ss_num(ss_i));
    dstn2  = sprintf([dstdir1 'KingKongWAM%05d.csv'], ss_num(ss_i));
    dstn21 = sprintf([dstdir1 'KingKongWAM%05d.bin'], ss_num(ss_i));
    dstn3  = sprintf([dstdir1 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn4  = sprintf([dstdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn5  = sprintf([dstdir1 'KingKongTSync.%05d.mat'], ss_num(ss_i));
%     dstn6  = sprintf([dstdir1 'KingKongEMG.%05d.csv'], ss_num(ss_i));
    dstn6  = sprintf([dstdir1 'KingKongEMG.%05d.mat'], ss_num(ss_i));
    dstn7  = sprintf([dstdir1 'KingKongOPT.%05d.csv'], ss_num(ss_i));
    
    dstnvd = sprintf([dstdirvd 'KingKong%05d.mp4'], ss_num(ss_i));
    dstnvdd= sprintf([dstdirvd 'KingKong%05d.txt'], ss_num(ss_i));
    
    try
        copyfile(fname1, dstn1);
    catch
    end
    
    try
        copyfile(fname2, dstn2);
    catch
    end
    
    try 
        copyfile(fname21, dstn21);
    catch
    end
    
    try
        copyfile(fname3, dstn3);
    catch
    end
    
    try
        copyfile(fname4, dstn4);
    catch
    end
    
    try 
        copyfile(fname5, dstn5);
    catch 
    end
        
    try 
        copyfile(fname6, dstn6);
    catch 
        try
        copyfile(fname60, dstn6);
        catch 
            try
                copyfile(fname61, dstn6);
            catch
%                 disp('no EMG at all');
            end
        end
    end
    
    try 
        copyfile(fname7, dstn7);
    catch 
        try
            copyfile(fname70, dstn7);
        catch 
%             disp('no OPT at all');
        end
    end
    
    % videologger
    try 
        copyfile(fnamevd, dstnvd);
        copyfile(fnamevdd, dstnvdd);
    catch 
%         disp('no Video in this session!');
    end

    disp(['session ' num2str(ss_num(ss_i)) ' copy condition:']);
    fprintf('||%s | %s |  %s | %s |  %s | %s | %s | %s ||\n', 'FMT', 'INT', 'FT', 'WAM', 'TS', 'EMG', 'OPT', 'VDL');
    fprintf('||%s   | %s   | %s   | %s   | %s   | %s   | %s   | %s   ||\n', cp_flag_str(1), cp_flag_str(2), cp_flag_str(3), cp_flag_str(4), cp_flag_str(5), cp_flag_str(6), cp_flag_str(7), cp_flag_str(8));
end


end

