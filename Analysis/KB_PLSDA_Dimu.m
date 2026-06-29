cd('C:\Users\Dimash\OneDrive\Рабочий стол\Thesis codebase.worktrees\First attempt worktree\Analysis');

clear vars
load('Merged_VNIR_SWIR_Summary.mat');

% Choosing the spectral range for analysis: VNIR, SWIR, or combined
X_vnir = cell2mat(T_merged.VNIR_Mean_Spectrum);
X_swir = cell2mat(T_merged.SWIR_Mean_Spectrum);
X_combined = cell2mat(T_merged.Combined_Spectrum);

% Choosing specific range 
X = X_combined; 

% Extracting the labels 
class_text = regexp(T_merged.Sample_ID, '^[^_]+', 'match', 'once');
class_cat = categorical(class_text);
class_names = categories(class_cat);
Y = double(class_cat);

%% Models for germination/no germination

rng(1); 

cal_idx = false(height(T_merged),1);
test_idx = false(height(T_merged),1);

for k = 1:numel(class_names)
    idx = find(Y == k);
    idx = idx(randperm(numel(idx))); 
    n_cal = round(0.7 * numel(idx));

    cal_idx(idx(1:n_cal)) = true; 
    test_idx(idx(n_cal+1:end)) = true;
end

c_data = X(cal_idx,:);
c_class = Y(cal_idx);
t_data = X(test_idx,:);
t_class = Y(test_idx);

[output_PLSDA, pretr, nlv] = PLSDA(c_data, c_class, t_data, t_class, 20, 100);

output_PLSDA{1,pretr}.SS{1,nlv}.Confusion_matrix
output_PLSDA{1,pretr}.SS{1,nlv}.CC
output_PLSDA{1,pretr}.SS{1,nlv}.Perclass_mean_accuracy



% model for old/new 
%% Models for old/new
% note unblanaced dataset: far more germinated than didn't
c_data = []; 
c_class= [];
c_sample=[];
for i = 1:size(cal,2)
x = cal(i).Mean_swir;
y = cal(i).age;
sam= cal(i).sample_no;
c_data = [c_data;x];
c_class = [c_class;y];
c_sample=[c_sample;sam]
end
clear x y

t_data = []; 
t_class= [];
t_sample=[];
for  i = 1:size(test,2)
x = test(i).Mean_swir;
y = test(i).age;
sam= test(i).sample_no;

t_data = [t_data;x];
t_class = [t_class;y];
t_sample=[t_sample;sam]
end
clear x y


[output_PLSDA_Age,pretr,nlv]=PLSDA(c_data,c_class,t_data,t_class,20,100);

%confusion matrix
output_PLSDA_Age{1,pretr}.SS{1,nlv}.Confusion_matrix
%overall accuracy
output_PLSDA_Age{1,pretr}.SS{1,nlv}.CC
%Mean class accuracy
output_PLSDA_Age{1,pretr}.SS{1,nlv}.Perclass_mean_accuracy