function [cp_flag] = cpDatarg2(ss_num)
%COPOYDATARG2 copy data from rg2 to my current computer with fixed source and desiny location 
%   srcDir: '/Volumes/rg2/data/KingKong/...'
%   dstDir: '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data'
%   copy may have one or multiple files. 

srcdir1 = '/Volumes/rg2/data/KingKong/Raw/';
srcdir2 = '/Volumes/rg2/data/KingKong/Formatted/';
srcdir3 = '/Volumes/rg2/data/KingKong/Intermediate/';

dstdir1 = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/';
dstdir2 = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/';

ss_num = ss_num(:);
cp_flag = zeros(size(ss_num));

ERRCODE_LOSTFORMAT = 8;
ERRCODE_LOSTIMTERM = 4;
ERRORCODE_LOSTFTCSV = 2;
ERRORCODE_LOSTWAMCSV= 1;
ERRCODE_LOSTTIMSYN = 16;


for ss_i = 1:length(ss_num)
    sprintf('copying... ');
    sprintf('%d / %d  ', ss_i, length(ss_num));
    
    % formatted data
    fname3 = sprintf([srcdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    if (~exist(fname3, 'file'))
        disp(['no formatted *.mat file for ss' num2str(ss_num(ss_i))]);
        cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTFORMAT;
    end
    
    % intermediate data
    fname4 = sprintf([srcdir3 'KingKong.%05d.mat'], ss_num(ss_i));
    if (~exist(fname4, 'file'))
        disp(['no intermediate *.mat file for ss' num2str(ss_num(ss_i))]);
        cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTIMTERM;
    end
    
    % FT data
    fname1 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname10= sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT0.csv'], ss_num(ss_i));
    if (~exist(fname1, 'file'))
        disp(['no processed FT file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname10, 'file'))
            disp('cautiion: no FT file at all');
            cp_flag(ss_i) = cp_flag(ss_i) + ERRORCODE_LOSTFTCSV;
        else
            disp('use FT0.csv file.');
        end
    end
    
    % WAM data
    fname2 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname20= sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM0.csv'], ss_num(ss_i));
    fname21= sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAMbin'], ss_num(ss_i));
    if (~exist(fname2, 'file'))
        disp(['no processed WAM file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname20, 'file'))
            disp('cautiion: no WAM file at all');
            cp_flag(ss_i) = cp_flag(ss_i) + ERRORCODE_LOSTWAMCSV;
        else
            disp('use WAM0.csv file.');
        end
        if (exist(fname21, 'file')) % binary file 
            disp('cautiion: WAM file is binary');
        end
    end
    
    % timeSync data
    fname5 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongTimeSync0.mat'], ss_num(ss_i));
    if (~exist(fname4, 'file'))
        disp(['no TimeSync *.mat file for ss' num2str(ss_num(ss_i))]);
        cp_flag(ss_i) = cp_flag(ss_i) + ERRCODE_LOSTTIMSYN;
    end
    
    % EMG data
    fname6 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongEMG%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname60= sprintf([srcdir1 'KingKong.DK.%05d/KingKongEMG0.csv'], ss_num(ss_i));
    if (~exist(fname6, 'file'))
        disp(['no processed EMG file for ss' num2str(ss_num(ss_i))]);
        if (~exist(fname60, 'file'))
            disp('cautiion: no EMG file at all');
%             cp_flag(ss_i) = cp_flag(ss_i) + ERRORCODE_LOSTWAMCSV;
        else
            disp('use EMG0.csv file.');
        end
    end
    
    dstn1  = sprintf([dstdir1 'KingKongFT%05d.csv'], ss_num(ss_i));
    dstn2  = sprintf([dstdir1 'KingKongWAM%05d.csv'], ss_num(ss_i));
    dstn21 = sprintf([dstdir1 'KingKongWAM%05d.bin'], ss_num(ss_i));
    dstn3  = sprintf([dstdir1 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn4  = sprintf([dstdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn5  = sprintf([dstdir1 'KingKongTSync.%05d.mat'], ss_num(ss_i));
    dstn6 = sprintf([dstdir1 'KingKongEMG.%05d.csv'], ss_num(ss_i));
    
    
    
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
        copyfile(fname60, dstn6);
    catch 
        
    end
end


end

