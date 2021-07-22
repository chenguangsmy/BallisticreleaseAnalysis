% Perturbation sequence generator 
% given 60 trials, make rule: 
%               5 - trials on turn;
%               12 - "turns" in total;
%               no perturbation at first 2 trials
%               eigher 2 or 3 trials catch trials in 5

% 1. times: variable describe: [if_pert, when_pert, when_release]
times = [
	0		1		3
	0		1		3
	0		1		3
	0		2		5
	0		3		5
	0		3		5
	0		3		7
	0		3		7
	0		3		7
	0		6		9
	1		1		3
	1		1		3
	1		1		3
	1		2		5
	1		3		5
	1		3		5
	1		3		7
	1		3		7
	1		3		7
	1		6		9];


% 2. combos: generate combos that satisfy the aforementioned criteria
    % number-of-random
    % which are 2 and which are 3, generate 6 1s out of 12 opsitions
    
    % generate positions
    combos = [];
    combo_num = 60;
    y_turn = sort(randsample(12, 6));
    for turn_i = 1:12
        combo_i = [];
        if turn_i == 1 % if turn(1), random 2/3
            if ismember(turn_i, y_turn) % all 3 perturb 
                
                combo_i = [(randsample(10,2)); (10 + randsample(10,3))];
            else % only 2 perturb
                combo_i = [randsample(10,2)];
                y = sort(randsample(3, 2));
                for i = 1:3
                    if ismember(i,y)
                        combo_i = [combo_i; (10 + randsample(10,1))]; 
                    else
                        combo_i = [combo_i; randsample(10,1)]; 
                    end
                end
            end
        else % other turn_i
            if ~ismember(turn_i, y_turn)     % else if (perturb_num == 2), random 2/5
                y = sort(randsample(5, 2));
                for i = 1:5
                    if ismember(i,y)
                        combo_i = [combo_i; (10 + randsample(10,1))]; 
                    else
                        combo_i = [combo_i; randsample(10,1)]; 
                    end
                end
            else % if (3), random 3/5
                y = sort(randsample(5, 3));
                for i = 1:5
                    if ismember(i,y)
                        combo_i = [combo_i; (10 + randsample(10,1))]; 
                    else
                        combo_i = [combo_i; randsample(10,1)]; 
                    end
                end
            end
        end
        combos = [combos combo_i'];
    end

perturb_schedue = reshape(combos>10, 5, 12)';
imagesc(perturb_schedue);

%combos 
if ~exist('targets_tmp.txt', 'file')
    fid = fopen( 'results.txt', 'wt' );
else
    fid = fopen( 'results.txt', 'wt' );
end

fprintf(fid, 'combos = [');
for combos_i = 1 : 60
    formatSpec = '%d  ';
    fprintf(fid, formatSpec, combos(combos_i));
    if mod(combos_i, 5) == 0
        fprintf(fid, '  ');
    end
end
fprintf(fid, '];');

fclose(fid);
