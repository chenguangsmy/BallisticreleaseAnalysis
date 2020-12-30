% Align trials 
%   Align by certain evident: STATES or time;
%   Align all trial elements (with index/time) for future use;
%   choose aligned array length (min/max/fixed)
function trial_ = ballisticRelease_aligntrials(trial, alignst, alignl)
% inputs: 
%   trial   the original trial structure
%   alignst the align state, can choice from ST_*
%   aligndl the align length, 'min', 'max', 'fix'

% Cleave 2020-02-08

ST_BGN = 1; 
ST_PRT = 2;
ST_FCR = 3;
ST_MOV = 4;
ST_HLD = 5;
ST_END = 6;
ST_RST = 7;

% original data
% e.g. align by ST_MOV
if nargin < 1,
    error('must have trial struct as input'),
elseif  (nargin < 2),
    alignst = ST_MOV;
    alignedl = 'max'; % aligned label
end

alignedfix = 1000; %

% for each trial
% find qualified trials 
qualified_trials = zeros(length(trial),1);
for ii = 1:length(trial),
    qualified_trials(ii) = sum(alignst == trial(ii).states);
end
trial = trial(find(qualified_trials));

% array length
switch alignl
    case 'min',
        alignedfix = min([trial.end]-[trial.bgn]);
    case 'max',
        alignedfix = max([trial.end]-[trial.bgn]);
    case 'fix',
        alignedfix = 1000;
end

% aligned-pt choose
if mod(alignedfix, 2) ~= 0,
    alignedfix = alignedfix + 1;
end
alignedpt = alignedfix/2 + 1; % index in trial_.* (post-align centre)

% aligned variables choose
fields = fieldnames(trial);
fields_idx = zeros(1,length(fields));
% search key words, which were contained in fields should be align.
keywords = {'time', 'position', 'force', 'Torque', 'vel'};
for i = keywords,
    cmp = strfind(fields, i{1});
    for subfield_i = 1:length(fields),
        fields_idx(subfield_i) = fields_idx(subfield_i) + ~isempty(cmp{subfield_i});
    end
end
fields_idx = find(fields_idx~=0);

trial_ = trial;
trial_num = length(trial);

% for each trial
for trial_i = 1:trial_num,
    event_index = [trial(trial_i).states_idx(trial(trial_i).states == alignst)]; % index in trial.* (pre-align centre)
    shift = alignedpt - event_index;
    % align state_index
    trial_(trial_i).states_idx = trial_(trial_i).states_idx + shift;
    % find range to avoid error
    event_index_bg = max(event_index-alignedfix/2, 1); % index no less than 1
    event_index_ed = min(event_index+alignedfix/2-1, length([trial(trial_i).time])); % no more than length
    event_index_l = event_index_bg - event_index; % datapoint before event_index in trial.*
    event_index_r = event_index_ed - event_index;
    
    % align for each variable
    for field_i = fields_idx,
        % dimension after aligned
        field_rows = size(trial(1).(fields{field_i}), 1);
        field_cols = alignedfix;
        
        % asign with nans
        eval(['trial_(trial_i).' fields{field_i} '= nan(field_rows, field_cols);']);

        % corresponds into trial_.*
        eval(['trial_(trial_i).' fields{field_i} '(:, alignedpt+event_index_l : alignedpt+event_index_r)' ...
            ' = trial(trial_i).' fields{field_i} '(:, event_index_bg : event_index_ed);']);
    end
end

end