% Currently I save the EMG using .txt file, where the data was written in
% string. This takes extra space in disk and process time. 
% I'm thinking if use .mat file or other file format will be a better
% solution. 
% cg, 2022-06-24 

%% 1. Example of the current data
clear; close all; clc;

sstmpemg = SessionScanEMG(4263);

t = sstmpemg.brtime;
dat = sstmpemg.dat;

tic
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_EMG/investigatingEMGsaver/dattmp4263.mat', 't', 'dat','-v6')
toc

% time: 0.029190 sec
% size: 8.9M
% dat_time: 368.89s (6min)
% dat_freq: 1000Hz

% Hence, the normal data (~320 trials, 2440s, 40min) will be ~60MB, and a 
%           save time of 0.19s
%        the most complicated data (~30*9*4=1080 trials, 7884s, 141min)
%           will be ~190MB, and a save time of 0.6198s

% totally fine saving time!!!

