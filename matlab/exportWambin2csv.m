function exportWambin2csv(name_bin, namecsv)
% function exportWambin2csv(name_bin, ss_num)
% this file is to deal with the WAM crash exception, convert saved binary
% file into *.csv file
% export the binary file to csv
% read Wam binary files and convert it into .csv file

% specified format as:
% |varname      | format        | size(byte) |
% | ----------- | ------------- | ---------- |
% |time         |     double    |     1*8    |
% |jointpos     | 4*1 double    |    DOF*8   |
% |jointvel     | 4*1 double    |    DOF*8   |
% |cartpos      | 3*1 double    |     3*8    |
% |cartvel      | 3*1 double    |     3*8    |
% |jointtqr     | 4*1 double    |    DO*8    |
% |cartfce      | 3*1 double    |     3*8    |
% |iteration    |       int     |     1*4    |
% |readtime     |       int     |     1*4    |
% |taskstate    |       int     |     1*4    |

% 0. initialize the conversion
DOF = 4;
data_ptn = {...% pattern
    % |varname  | format       | #length |  #byte |
    'time',     'double',       1,          8,
    'jp',       'double',       DOF,        8,
    'jv',       'double',       DOF,        8,
    'cp',       'double',       3,          8,
    'cv',       'double',       3,          8,
    'jt',       'double',       DOF,        8,
    'cf',       'double',       3,          8,
    'it',       'int',          1,          4,
    'rdt',      'int',          1,          4,
    'ts'        'int',          1,          4,
    };

datwidth = sum([data_ptn{:,3}].*[data_ptn{:,4}]);
% specify the file name
fname = name_bin; %'bt202009046AuSVW';
fnamew= namecsv;%'exportWAM1.csv';

% 1. form a mat saving variables
fileID = fopen(fname,'r');
datatmp = fread(fileID);
datalen = ceil(size(datatmp,1)/datwidth);
datamat = zeros(datalen,sum([data_ptn{:,3}]));
clear datatmp

% 2. read varible from binary file
fileID = fopen(fname,'r');
n = 0;
while ~feof(fileID)
    n = n+1;
    idx = 1;
    for ri = 1:size(data_ptn,1)
        idx_stt = idx;
        idx_edn = idx + data_ptn{ri,3}-1;
        data = fread(fileID,data_ptn{ri,3},data_ptn{ri,2});
        if (n>datalen-1)
            continue % read out
        end
        datamat(n,idx:idx_edn) = data;
        idx = idx_edn + 1;
        %eval([data_ptn{ri,1} '=transpose(data);']);
    end
end
fclose(fileID);

% 3. export variable to csv file
write_spec = [];
for ri = 1:size(data_ptn, 1)
    if strcmp(data_ptn{ri,2}, 'double')
        write_type = '%f';
    elseif strcmp(data_ptn{ri,2}, 'int')
        write_type = '%d';
    end
    sep_type = ','; % csv
    if ri == size(data_ptn,1)
        sep_type = '\n';
    end
    %data = eval(data_ptn{ri,1});
    for di = 1:data_ptn{ri,3}
        write_spec = [write_spec, write_type, sep_type];
    end
end

wfileID = fopen(fnamew,'w');
for ti = 1:size(datamat,1)
    fprintf(wfileID, write_spec, datamat(ti,:));
end

fclose(wfileID);

end