
%% Plotting the validation accuracy and the final confusion matrix.
function [output_PLSDA, opt_pretr, opt_nlv, Fold_Mean_Class_Accuracy] = ...
    PLSDA_new(Xcal, Ycal, Xtest, Ytest, maxnlv, niters)

    % Choosing the preprocessing path and LV count from the calibration folds.
    [~, ~, Fold_Mean_Class_Accuracy] = ...
        Npls_da_N_SEP_pre_randCV_SS(Xcal, Ycal, maxnlv, niters);

    [opt_pretr, opt_nlv] = selectOptimalPretreatmentAndLV(Fold_Mean_Class_Accuracy);

    % Fitting the final model on the full calibration set and scoring the hold-out test set.
    output_PLSDA = Npls_da_N_SEP_pre_SS(Xcal, Ycal, Xtest, Ytest, maxnlv);
    Pred_Ytest = output_PLSDA{opt_pretr}.Pred_class_t(:, opt_nlv);

    % Displaying the cross-validated accuracy curves and the test-set confusion matrix.
    figure('Name', 'PLS-DA model selection', 'Color', 'w');
    tiledlayout(1, 2);

    nexttile;
    plot(Fold_Mean_Class_Accuracy', '-o', 'LineWidth', 2);
    xlabel('Number of latent variables', 'FontSize', 20);
    ylabel('Mean class accuracy', 'FontSize', 20);
    legend({'none', 'SNV', '1D', '2D', 'SNV + 1D', 'SNV + 2D', 'MSC'}, ...
        'Location', 'best');
    ylim([0, 100]);
    ax = gca;
    ax.FontSize = 16;

    nexttile;
    confusionchart(double(output_PLSDA{1}.Yval), Pred_Ytest);
    title('Test set');
end

function [opt_pretr, opt_nlv] = selectOptimalPretreatmentAndLV(Fold_Mean_Class_Accuracy)
    nPretreatments = size(Fold_Mean_Class_Accuracy, 1);
    nLv = size(Fold_Mean_Class_Accuracy, 2);

    % Stopping when the accuracy gain drops below one percentage point.
    accuracyDifference = diff(Fold_Mean_Class_Accuracy, 1, 2);
    belowThreshold = accuracyDifference < 1;
    nLvLimit = zeros(nPretreatments, 1);

    for pretreatmentIndex = 1:nPretreatments
        firstBelowThreshold = find(belowThreshold(pretreatmentIndex, :), 1);
        if isempty(firstBelowThreshold)
            nLvLimit(pretreatmentIndex) = nLv;
        else
            nLvLimit(pretreatmentIndex) = firstBelowThreshold + 1;
        end
    end

    accuracyAtLimit = zeros(nPretreatments, 1);
    for pretreatmentIndex = 1:nPretreatments
        accuracyAtLimit(pretreatmentIndex) = ...
            Fold_Mean_Class_Accuracy(pretreatmentIndex, nLvLimit(pretreatmentIndex));
    end

    [~, opt_pretr] = max(accuracyAtLimit);
    opt_nlv = nLvLimit(opt_pretr);

    
    Fold_sub = Fold_Mean_Class_Accuracy(:, 1:opt_nlv);
    [pretrMax, nlvMax] = max(Fold_sub, [], 2);
    [~, opt_pretr] = max(pretrMax);
    opt_nlv = nlvMax(opt_pretr);
end


% Comparing the seven preprocessing paths across repeated random validation folds.
function [Fold_CV, Fold_SeSp, Fold_Mean_Class_Accuracy, Mean_Class_Accuracy, B_all] = ...
    Npls_da_N_SEP_pre_randCV_SS(Xcal, Ycal, maxnlv, niters)

    Xcal = double(Xcal);
    nPretreatments = 7;
    nSamples = size(Xcal, 1);

    CC_CV = zeros(nPretreatments, maxnlv, niters);
    SeSp_CV = zeros(nPretreatments, maxnlv, niters);
    Mean_Class_Accuracy = zeros(nPretreatments, maxnlv, niters);
    B_all = cell(niters, nPretreatments);

    for iteration = 1:niters
        rng(iteration);
        randomIndices = randperm(nSamples);
        validationCount = floor(nSamples * 0.3);

        Xcal2 = Xcal(randomIndices(1:validationCount), :);
        Ycal2 = Ycal(randomIndices(1:validationCount), :);
        Xcal1 = Xcal(randomIndices(validationCount + 1:end), :);
        Ycal1 = Ycal(randomIndices(validationCount + 1:end), :);

        output1 = Npls_da_N_SEP_pre_SS(Xcal1, Ycal1, Xcal2, Ycal2, maxnlv);

        for pretreatmentIndex = 1:nPretreatments
            CC_CV(pretreatmentIndex, :, iteration) = output1{pretreatmentIndex}.CC_test;
            SeSp_CV(pretreatmentIndex, :, iteration) = output1{pretreatmentIndex}.SeSp;

            for lvIndex = 1:maxnlv
                Mean_Class_Accuracy(pretreatmentIndex, lvIndex, iteration) = ...
                    output1{pretreatmentIndex}.SS{lvIndex}.Perclass_mean_accuracy;
            end

            B_all{iteration, pretreatmentIndex} = output1{pretreatmentIndex}.B;
        end
    end

    Fold_CV = mean(CC_CV, 3);
    Fold_SeSp = mean(SeSp_CV, 3);
    Fold_Mean_Class_Accuracy = mean(Mean_Class_Accuracy, 3);
end

% Applying the preprocessing variants before fitting the one-vs-rest PLS-DA model.
function output = Npls_da_N_SEP_pre_SS(Xcal, Ycal, Xval, Yval, maxnlv)
    Xcal = double(Xcal);
    Xval = double(Xval);

    output = cell(1, 7);
    output{1} = Npls_da_N_SEP_SS(Xcal, Ycal, Xval, Yval, maxnlv);
    output{2} = Npls_da_N_SEP_SS(SNV(Xcal), Ycal, SNV(Xval), Yval, maxnlv);
    output{3} = Npls_da_N_SEP_SS(savgol(Xcal, 15, 3, 1), Ycal, savgol(Xval, 15, 3, 1), Yval, maxnlv);
    output{4} = Npls_da_N_SEP_SS(savgol(Xcal, 15, 3, 2), Ycal, savgol(Xval, 15, 3, 2), Yval, maxnlv);
    output{5} = Npls_da_N_SEP_SS(savgol(SNV(Xcal), 15, 3, 1), Ycal, savgol(SNV(Xval), 15, 3, 1), Yval, maxnlv);
    output{6} = Npls_da_N_SEP_SS(savgol(SNV(Xcal), 15, 3, 2), Ycal, savgol(SNV(Xval), 15, 3, 2), Yval, maxnlv);
    output{7} = Npls_da_N_SEP_SS(MSC(Xcal, mean(Xcal)'), Ycal, MSC(Xval, mean(Xcal)'), Yval, maxnlv);
end
%% functions


function [output] = Npls_da_N_SEP_SS(Xcal,Ycal,Xval,Yval,maxnlv)
% Encoding each class in a one-vs-rest format and fitting a latent-variable classifier.
Ycal_unique = unique(Ycal);




nclasses = size(Ycal_unique, 1);
nCalibrationSamples = size(Ycal, 1);
nValidationSamples = size(Yval, 1);

Ycal_num = zeros(nCalibrationSamples, 1);
Yval_num = zeros(nValidationSamples, 1);
Ycal_cl = false(nCalibrationSamples, nclasses);
Yval_cl = false(nValidationSamples, nclasses);


for classIndex = 1:nclasses
    Ycal_cl(:, classIndex) = Ycal == Ycal_unique(classIndex);
    Yval_cl(:, classIndex) = Yval == Ycal_unique(classIndex);

    Ycal_num(Ycal == Ycal_unique(classIndex)) = classIndex;
    Yval_num(Yval == Ycal_unique(classIndex)) = classIndex;
end

nLatentVariables = maxnlv;
B = NplsN(Xcal, Ycal_cl, nLatentVariables);
Y_pred_cal_all = zeros(nCalibrationSamples, nLatentVariables, nclasses);
Y_pred_test_all = zeros(nValidationSamples, nLatentVariables, nclasses);
c0_m = zeros(nLatentVariables, nclasses);
c0_sd = zeros(nLatentVariables, nclasses);
c1_m = zeros(nLatentVariables, nclasses);
c1_sd = zeros(nLatentVariables, nclasses);
PO_Y = zeros(nCalibrationSamples, nLatentVariables, nclasses);
P1_Y = zeros(nCalibrationSamples, nLatentVariables, nclasses);
threshold = zeros(nLatentVariables, nclasses);
Err_per_cal = zeros(nLatentVariables, 1);
Err_per_test = zeros(nLatentVariables, 1);
P0 = zeros(1, nclasses);
P1 = zeros(1, nclasses);

for lvIndex = 1:nLatentVariables

Y_pred_cal = Xcal * B(:, :, lvIndex);
Y_pred_test = Xval * B(:, :, lvIndex);

% Estimating the Gaussian class margins and posterior probabilities for each LV step.

    Err_cal = false(nCalibrationSamples, nclasses);
    Err_test = false(nValidationSamples, nclasses);
    for classIndex = 1:nclasses

        Ind1_c = find(Ycal_cl(:, classIndex));
        Ind0_c = find(~Ycal_cl(:, classIndex));

        Y_pred_class1 = Y_pred_cal(Ind1_c, classIndex);
        Y_pred_class0 = Y_pred_cal(Ind0_c, classIndex);

        c0_m(lvIndex, classIndex) = mean(Y_pred_class0);
        c0_sd(lvIndex, classIndex) = std(Y_pred_class0);

        c1_m(lvIndex, classIndex) = mean(Y_pred_class1);
        c1_sd(lvIndex, classIndex) = std(Y_pred_class1);

% Calculating the Gaussian likelihoods of the positive and negative class scores.
        Py_0 = 1 ./ (sqrt(2*pi) * c0_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_cal(:, classIndex) - c0_m(lvIndex, classIndex)) ...
            / c0_sd(lvIndex, classIndex)).^2);
        Py_1 = 1 ./ (sqrt(2*pi) * c1_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_cal(:, classIndex) - c1_m(lvIndex, classIndex)) ...
            / c1_sd(lvIndex, classIndex)).^2);

% Weighting the class priors from the calibration split.
        P0(classIndex) = numel(Ind0_c) / (numel(Ind1_c) + numel(Ind0_c));
        P1(classIndex) = numel(Ind1_c) / (numel(Ind1_c) + numel(Ind0_c));

        P0_y = Py_0 * P0(classIndex) ./ ...
            (Py_0 * P0(classIndex) + Py_1 * P1(classIndex));
        P1_y = Py_1 * P1(classIndex) ./ ...
            (Py_0 * P0(classIndex) + Py_1 * P1(classIndex));


        PO_Y(:, lvIndex, classIndex) = P0_y;
        P1_Y(:, lvIndex, classIndex) = P1_y;
        Y_pred_cal_all(:, lvIndex, classIndex) = Y_pred_cal(:, classIndex);
        Y_pred_test_all(:, lvIndex, classIndex) = Y_pred_test(:, classIndex);

% Setting the decision boundary from the two Gaussian class distributions.

        threshold(lvIndex, classIndex) = ...
            (c1_m(lvIndex, classIndex) - 2*c1_sd(lvIndex, classIndex) ...
            - (c0_m(lvIndex, classIndex) + 2*c0_sd(lvIndex, classIndex))) / 2 ...
            + c0_m(lvIndex, classIndex) + 2*c0_sd(lvIndex, classIndex);

        Err_cal(:, classIndex) = ...
            (Y_pred_cal(:, classIndex) > threshold(lvIndex, classIndex)) ~= ...
            Ycal_cl(:, classIndex);
        Err_test(:, classIndex) = ...
            (Y_pred_test(:, classIndex) > threshold(lvIndex, classIndex)) ~= ...
            Yval_cl(:, classIndex);


    end
    Err_per_cal(lvIndex) = 100 * nnz(Err_cal) / numel(Err_cal);
    Err_per_test(lvIndex) = 100 * nnz(Err_test) / numel(Err_test);
end

%figure,plot(Err_per_cal,'o'),xlabel('Number of LV'),ylabel('% Error')

%maxnlv = input('Maximum number of LVs?: ', 's');
%maxnlv=str2double(maxnlv);
%maxnlv=20
P1_y_c = zeros(nCalibrationSamples, nclasses, nLatentVariables);
P1_y_t = zeros(nValidationSamples, nclasses, nLatentVariables);

for lvIndex = 1:nLatentVariables
    for classIndex = 1:nclasses

% Estimating the posterior probability for each class and selecting the most likely label.
        Py_0_c = 1 ./ (sqrt(2*pi) * c0_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_cal_all(:, lvIndex, classIndex) ...
            - c0_m(lvIndex, classIndex)) / c0_sd(lvIndex, classIndex)).^2);
        Py_1_c = 1 ./ (sqrt(2*pi) * c1_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_cal_all(:, lvIndex, classIndex) ...
            - c1_m(lvIndex, classIndex)) / c1_sd(lvIndex, classIndex)).^2);

        P1_y_c(:, classIndex, lvIndex) = Py_1_c * P1(classIndex) ./ ...
            (Py_0_c * P0(classIndex) + Py_1_c * P1(classIndex));

        Py_0_t = 1 ./ (sqrt(2*pi) * c0_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_test_all(:, lvIndex, classIndex) ...
            - c0_m(lvIndex, classIndex)) / c0_sd(lvIndex, classIndex)).^2);
        Py_1_t = 1 ./ (sqrt(2*pi) * c1_sd(lvIndex, classIndex)) * ...
            exp(-0.5 * ((Y_pred_test_all(:, lvIndex, classIndex) ...
            - c1_m(lvIndex, classIndex)) / c1_sd(lvIndex, classIndex)).^2);

        P1_y_t(:, classIndex, lvIndex) = Py_1_t * P1(classIndex) ./ ...
            (Py_0_t * P0(classIndex) + Py_1_t * P1(classIndex));

    end
end

[~,Pred_class_c] = max(P1_y_c, [], 2);
[~,Pred_class_t] = max(P1_y_t, [], 2);

cl_error_cal = squeeze(sum(Pred_class_c ~= Ycal_num) / nCalibrationSamples);
cl_error_test = squeeze(sum(Pred_class_t ~= Yval_num) / nValidationSamples);



CC_all_c = 100 .* (1 - cl_error_cal);
CC_all_t = 100 .* (1 - cl_error_test);

% Mapping the index-based predictions back to the original class labels.
Pred_class_c_actual = zeros(nCalibrationSamples, nLatentVariables, 'like', Ycal_unique);
Pred_class_t_actual = zeros(nValidationSamples, nLatentVariables, 'like', Ycal_unique);
for classIndex = 1:nclasses
    for lvIndex = 1:nLatentVariables
        Pred_class_c_actual(Pred_class_c(:, :, lvIndex) == classIndex, lvIndex) = ...
            Ycal_unique(classIndex);
        Pred_class_t_actual(Pred_class_t(:, :, lvIndex) == classIndex, lvIndex) = ...
            Ycal_unique(classIndex);
    end
end

% Summarizing sensitivity and specificity for each LV setting.
nValidationClasses = numel(unique(Yval));
SeSp = zeros(1, nLatentVariables);
SS = cell(1, nLatentVariables);
if nValidationClasses > 1
    for lvIndex = 1:nLatentVariables
        SS{lvIndex} = SensSpecConf_Nclass(Yval, Pred_class_t_actual(:, lvIndex));
        classSeSp = SS{lvIndex}.Sensitivity .* SS{lvIndex}.Specificity;
        SeSp(lvIndex) = mean(classSeSp);
    end
else
    for lvIndex = 1:nLatentVariables
        SS{lvIndex} = [Yval'; Pred_class_t_actual(:, lvIndex)'];
    end
end
output.SeSp = SeSp;
output.SS=SS;



% plot(100.*[cl_error_cal,cl_error_test],'o-')
% legend({'cal','test'})
% ylabel('Classification Error (%)')
% xlabel('Number of Latent Variables')
output.B=B;
output.PO_Y=PO_Y;
output.P1_Y=P1_Y;
output.Y_pred_cal=Y_pred_cal_all;
output.Y_pred_test=Y_pred_test_all;
output.c0_m=c0_m;
output.c0_sd=c0_sd;
output.c1_m=c1_m;
output.c1_sd=c1_sd;
output.threshold=threshold;
output.nclasses=nclasses;
output.P0=P0;
output.P1=P1;
output.nlv = nLatentVariables;
output.CC_test=CC_all_t;
output.CC_cal=CC_all_c;
output.Pred_class_c=Pred_class_c_actual;
output.Pred_class_t=Pred_class_t_actual;
output.Ycal=Ycal;
output.Yval=Yval;

end

function [B,T,TT,P,W,Q] = NplsN(X,Y,A)
% Updating the latent scores and regression weights with the NIPALS routine.
X = double(X);
Y = double(Y);
nPredictors = size(X, 2);
nResponses = size(Y, 2);
W = zeros(nPredictors, A);
T = zeros(size(X, 1), A);
P = zeros(nPredictors, A);
Q = zeros(A, nResponses);
B = zeros(nPredictors, nResponses, A);

for lvIndex = 1:A
    v = X' * Y(:, 1);
    W(:, lvIndex) = v / sqrt(v' * v);
    T(:, lvIndex) = X * W(:, lvIndex);
    TT = T(:, lvIndex)' * T(:, lvIndex);
    P(:, lvIndex) = X' * T(:, lvIndex) / TT;
    X = X - T(:, lvIndex) * P(:, lvIndex)';
    Q(lvIndex, :) = T(:, lvIndex)' * Y / TT;
    B(:, :, lvIndex) = W / (P' * W) * Q;
end
end

function[output]=SensSpecConf_Nclass(Actual,Predicted)
% Building the confusion matrix and summarizing class-level metrics.

nclasses = unique(Actual);
nClassCount = numel(nclasses);
Confusion_matrix = zeros(nClassCount, nClassCount);

for actualIndex = 1:nClassCount
    actualMask = Actual == nclasses(actualIndex);
    for predictedIndex = 1:nClassCount
        Confusion_matrix(actualIndex, predictedIndex) = ...
            sum(Predicted(actualMask) == nclasses(predictedIndex));
    end
end

CC = 100 * sum(diag(Confusion_matrix)) / sum(Confusion_matrix(:));

% Calculating the per-class accuracy from the confusion matrix.
Per_class_accuracy = zeros(1, nClassCount);
for classIndex = 1:nClassCount
    actualMask = Actual == nclasses(classIndex);
    Per_class_accuracy(classIndex) = ...
        100 * sum(Predicted(actualMask) == nclasses(classIndex)) / nnz(actualMask);
end

% Computing sensitivity, specificity, precision, and F1 across the classes.
TP_class = diag(Confusion_matrix)';
FN_class = sum(Confusion_matrix, 2)' - TP_class;
FP_class = sum(Confusion_matrix, 1) - TP_class;
TN_class = sum(Confusion_matrix(:)) - TP_class - FN_class - FP_class;

Sensitivity = TP_class ./ (TP_class + FN_class);

Specificity = TN_class ./ (FP_class + TN_class);

Precision = TP_class ./ (TP_class + FP_class);

F1 = 2 * TP_class ./ (2 * TP_class + FP_class + FN_class);

output.Confusion_matrix=Confusion_matrix;
output.CC=CC;
output.Sensitivity=Sensitivity;
output.Specificity=Specificity;
output.Precision=Precision;
output.F1 = F1;
output.F1mean=mean(F1);
output.Perclass_mean_accuracy=mean(Per_class_accuracy);      
end

function Xsnv=SNV(X)
% Centering and scaling each spectrum to unit variance.
X = double(X);
Xsnv = (X - mean(X,2))./(std(X,0,2));
Xsnv(isnan(Xsnv)) = 0;
end

function [z]=savgol(x,F,K,Dn)
% Applying a Savitzky-Golay filter and derivative across each spectrum row.

if nargin < 4 || isempty(Dn)
    Dn = 0;
end

x = double(x);
aux=F;
gap=floor(F/2);
if(gap==F/2)
    F=F+1;
    warning('savgol:EvenWindowAdjusted', 'Window size must be odd. Using F = %d instead of %d.', F, aux);
end

[~,g]=sgolaycoef(K,F);
[nrow,ncol]=size(x);
M=zeros(nrow,ncol);

halfWindow = (F-1)/2;
leftPad = repmat(x(:,1),1,halfWindow);
rightPad = repmat(x(:,end),1,halfWindow);
xPad = [leftPad, x, rightPad];

for i=1:nrow
    y=xPad(i,:);
    for j=1:ncol
        idx=j:j+F-1;
        if Dn==0
            M(i,j)=g(:,1)'*y(idx)';
        else
            M(i,j)=Dn*g(:,Dn+1)'*y(idx)';
        end
    end
end

z=M;
end

%%%%%%%%%%%%%
function [B,G] = sgolaycoef(k,F)
% Generating the polynomial coefficients used for smoothing and differentiation.

W = eye(F);
halfWindow = (F - 1) / 2;
s = fliplr(vander(-halfWindow:halfWindow));
S = s(:,1:k+1);

[~,R] = qr(sqrt(W) * S,0);

G = (R \ (R' \ S'))';

B = G*S'*W;

end

function[Xcor]=MSC(X,Target_mean)
% Removing multiplicative scatter effects by fitting a linear regression correction.
Xcal=[ones(size(Target_mean)),Target_mean];
Xcor=zeros(size(X));
for i=1:size(X,1)
Ycal=X(i,:)';   
b = regress(Ycal,Xcal);
Xcor(i,:)=((Ycal-b(1)))/b(2);
end
end