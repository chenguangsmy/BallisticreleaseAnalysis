%% 
% get the maximum voluntary force from these blocks.. 
% MVF will be shown as the experiment sequence: right, left, front, back
% 
% |subject| LR    | FB     |
% |------| ------ | ------ |
% |subj1 | 4399   | 4403   |
% |subj2 | 4407   | 4411   |
% |subj3 | 4417   | 4420   |
% |subj4 | 4431   | 4435   |
% |subj5 |        |        | 4284
% |subj6 | 4311   | 4315   |
% |subj7 | 4445   | 4449   |
% |subj8 | 4454   | 4457   |
% |subj9 | 4462   | 4465   |
% |subj10| 4471   | 4475   |
% |subj11| 4480   | 4484   |
% |subj12| 4490   | 4493   |
% |subj13| 4499   | 4502   |
% |subj14| 4511   | 4514   |
% |subj15| 4519   | 4522   |
% |subj16| 4529   | 4532   |
% |subj17| 4540   | 4544   |
% |subj18| 4556   | 4561   |
% |subj19| 4567   | ----   |
% |subj20| 4572   | 4575   |

% the force sequence should be: right-left-front-back
% the data should be a schalar in each direction...

ss_num = [%  4399    4403   
            4407    4411   
            4417    4420   
            4431    4435   
            4284    4284   
            4315    4312   
            4445    4449   
            4454    4457   
            4462    4465   
            4470    4475   
            4480    4484   
            4490    4493   
            4499    4502   
            4511    4514   
            4519    4522   
            4529    4532   
            4540    4544   
            4556    4561   
            4567    nan   
            4572    4575];

% % % check the time of take average 
% % ss_num_line = ss_num(1:end,:)';
% % ss_num_line = ss_num_line(:);
% % ss_num_line = ss_num_line([1:end-3, end-1:end])
% % 
% % for i = 1:length(ss_num_line)
% % figure(i)
% % i
% % fttmp = SessionScanFT(ss_num_line(i));
% % fttmp.plotForceOrigin(); 
% % title(['session' num2str(ss_num_line(i))]);
% % i = i + 1;
% % end

time_slots = {% 1st subject
    [63580       69190       76240       78860      185280      195410      205140      214630] [ 78944       87754       95945      103848      112154      120193      125140      132802]
    [17029       25141       38266       45810       61067       69069       78961       85958] [16208       24343       38138       45846       62430       70689       85431       95094]
    [220244      227505      252413      260591      277182      287718      304429      313359] [313930      322172      337351      346114      362306      371289      386734      395761]
    [19379       25230       32840       38246       45671       51801       59372       65445] [123394      128706      139962      146368      157828      166415      178417      184501]
    [445279      453871      463028      473634      496643      505494      512928      521589] [236949      246703      267859      275727      292557      299800      312827      321104] 
    [287527      297843      330894      340597      368118      377403      418477      425805] [164382      174540      200781      209887      240041      248529      266527      275044]
    [181196      190016      208304      216717      233155      241095      259709      268235] [360771      369201      388682      396644      414724      423444      437852      448466]
    [203028      212347      227830      237899      253898      262282      280114      291283] [188413      198082      219649      230019      246307      256220      271208      280225]
    [129646      139555      159496      168065      188362      197769      213075      221324] [23284       31504       47696       56566       79841       89170      106221      113997]
    [261615      272201      290237      300602      326007      335173      355171      365388] [127764      136115      153324      162320      179788      188286      204332      211250]
    [297371      307251      325095      333980      351062      360200      375714      385342] [9987       18562       34030       42370       62553       69654       84445       92343]
    [328566      340331      355527      364644      385328      395336      409741      418964] [108687      117833      132400      141816      157324      166308      180769      190392]
    [119038      128921      144117      154454      170205      180464      199083      207788] [143632      152361      167870      177526      193034      201618      227541      237704]
    [178736      187386      207965      217537      243496      253139      270822      280111] [18017       26699       42802       52388       68693       78884       93277      101976]
    [143590      152512      168785      176586      201007      209363      233894      242832] [79756       89099      105462      114200      129534      139112      153167      162371]
    [221366      231290      246076      254493      269356      279972      292959      302887] [77908       87795      105630      114935      135080      143806      168635      176627]
    [360839      369980      394009      405575      422586      433423      452340      461662] [100099      108981      100099      108981      164442      174023      198606      208701]
    [91897      102442      118724      127984      146172      155218      170217      179210] [nan nan nan nan nan nan nan nan]
    [252914      264722      283011      293735      311890      321907      342726      354081] [70894       81298       98516      108586      137849      148119      165112      174342]
    }


% t_all = [252914	-57.112705
% 264722	-51.180939
% 283011	-51.084083
% 293735	-45.286692
% 311890	74.499811
% 321907	64.099693
% 342726	71.390665
% 354081	54.474843
% ];
% t_all(:,1)'

%% for each time slot pair, do the alculating of average... 
% 
mvf_mat = zeros(19, 4); % for each direction, the mvf is one value
for subj_i = 1:19 
    rot_i = 1; % left-right 
    fttmp = SessionScanFT(ss_num(subj_i,rot_i));
    fce_0 = mean(fttmp.force_origin(:,end-100:end), 2); 
    % get each individual force...
    fce_idx11 = time_slots{subj_i,1}(1):time_slots{subj_i,1}(2);
    fce11 = vecnorm(fttmp.force_origin(:,fce_idx11) - fce_0);
    fce11 = mean(fce11); 

    fce_idx12 = time_slots{subj_i,1}(3):time_slots{subj_i,1}(4);
    fce12 = vecnorm(fttmp.force_origin(:,fce_idx12) - fce_0);
    fce12 = mean(fce12); 

    fce_idx21 = time_slots{subj_i,1}(5):time_slots{subj_i,1}(6);
    fce21 = vecnorm(fttmp.force_origin(:,fce_idx21) - fce_0);
    fce21 = mean(fce21); 

    fce_idx22 = time_slots{subj_i,1}(7):time_slots{subj_i,1}(8);
    fce22 = vecnorm(fttmp.force_origin(:,fce_idx22) - fce_0);
    fce22 = mean(fce22); 
    mvf_mat(subj_i,1) = (fce11+fce12)/2;
    mvf_mat(subj_i,2) = (fce21+fce22)/2;

    figure(); hold on; 
    plot([fttmp.force_origin - fce_0]');
    f_avg = mean(fttmp.force_origin(:,fce_idx11) - fce_0,2);
    yline(f_avg, 'linewidth', 2, 'Color',[0.5 0.5 0.5]);
    title(['ss\_' num2str(ss_num(subj_i,1))]);
        
    rot_i = 2; % front-back 
    if subj_i == 18
        continue
    end
    fttmp = SessionScanFT(ss_num(subj_i,rot_i));
    fce_idx11 = time_slots{subj_i,2}(1):time_slots{subj_i,2}(2);
    fce11 = vecnorm(fttmp.force_origin(:,fce_idx11) - fce_0);
    fce11 = mean(fce11); 

    fce_idx12 = time_slots{subj_i,2}(3):time_slots{subj_i,2}(4);
    fce12 = vecnorm(fttmp.force_origin(:,fce_idx12) - fce_0);
    fce12 = mean(fce12); 

    fce_idx21 = time_slots{subj_i,2}(5):time_slots{subj_i,2}(6);
    fce21 = vecnorm(fttmp.force_origin(:,fce_idx21) - fce_0);
    fce21 = mean(fce21); 

    fce_idx22 = time_slots{subj_i,2}(7):time_slots{subj_i,2}(8);
    fce22 = vecnorm(fttmp.force_origin(:,fce_idx22) - fce_0);
    fce22 = mean(fce22); 

    mvf_mat(subj_i,3) = (fce11+fce12)/2;
    mvf_mat(subj_i,4) = (fce21+fce22)/2;

    figure(); hold on; 
    plot([fttmp.force_origin - fce_0]');
    f_avg = mean(fttmp.force_origin(:,fce_idx11) - fce_0,2);
    yline(f_avg, 'linewidth', 2, 'Color',[0.5 0.5 0.5]);
    title(['ss\_' num2str(ss_num(subj_i,1))]);

end

subj_idx = [1:20];
mvf_mat = [nan nan nan nan; mvf_mat];
mvf_mat(19,3:4) = nan;
save('mvf_matfile.mat', 'mvf_mat', 'subj_idx');

% plot the mvf in the bar plot
subj_num = 20;
figure('name', 'MVF'); hold on;
bar(1:4, mean(mvf_mat(:,:),1, 'omitnan'));
errorbar(mean(mvf_mat(:,:),1, 'omitnan'), std(mvf_mat(:,:),1, 'omitnan'), 'LineWidth',3, 'Color',[0 0 0]); 
x_randomnized_mat = [1:4]+(rand(subj_num,4)-0.5)/2;
mvf_all_sel = mvf_mat(subj_idx,:);
x_randomnized_mat_sel = x_randomnized_mat(subj_idx,:);
scatter(x_randomnized_mat_sel(:), mvf_all_sel(:), 'b'); 
xticks([1:4]);
% ylim([0 1]);
% xticklabels({'15N 2.5cm', '15N 5.0cm', '15N 7.5cm', '20N 2.5cm', '20N 5.0cm', '20N 7.5cm', '25N 2.5cm', '25N 5.0cm', '25N 7.5cm'});
xticklabels({'1 (right)', '2 (left)', '3 (front)', '4 (back)'});
xlabel('directions'); 
ylabel('force (N)');
legend('averaged force','force std','single subject');
title({'Maximum Voluntary Force across directions'}); 

%% check the mvf relationship with the body mass
mass_all = [63.2000   65.7000   85.8000   50.2000   82.0000   77.8000   90.0000   78.5000   81.3000   60.0000   75.0000  113.0000  101.5000   55.2000   65.0000   62.1000  58.0000   58.5000   52.2000   57.3000];
mass_mat = repmat(mass_all, 4,1)';

figure('name', 'MVF vs mass');
scatter(mass_mat(:), mvf_mat(:), 'black');
xlabel('mass (kg)'); 
ylabel('mvf (N)');
refline
title('mass with MVF');

% seperate male and female - not a good figure, abort!
gender_all = 'MMMMMMMMMFMMMFMFFFFF'; 
male_idx = gender_all == 'M';
female_idx = gender_all == 'F';

mass_male = mass_mat(male_idx,:);
mvf_male =  mvf_mat(male_idx,:);
mass_female = mass_mat(female_idx,:);
mvf_female =  mvf_mat(female_idx,:);
figure('name', 'MVF vs mass'); hold on; 
scatter(mass_male(:), mvf_male(:), 'blue');
scatter(mass_female(:), mvf_female(:), 'red');
xlabel('mass (kg)'); 
ylabel('mvf (N)');