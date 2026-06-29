%% AG June 2026 
% PLSDA with pretreatments, given a Cal and Test set, do random Cross-Val, optimise and predict on the test set
%% plot output as Cross validation Mean Class Accuracy and Confusion matrix of test set
function [output_PLSDA,opt_pretr,opt_nlv]=PLSDA(Xcal,Ycal,Xtest,Ytest,maxnlv,niters)
[Fold_CV, Fold_SeSp, Fold_Mean_Class_Accuracy,Mean_Class_Accuracy,B_all]=Npls_da_N_SEP_pre_randCV_SS(Xcal,Ycal,maxnlv,niters)

 %calculate difference in accuracy for different nLV
       Acc_diff=Fold_Mean_Class_Accuracy(:,3:20)-Fold_Mean_Class_Accuracy(:,2:19);
       %select nLV for each pretr for which diff < 1
       [~,inLV]=max((Acc_diff<1)');
       for ipretr=1:7
           Acc_nLV(ipretr)=Fold_Mean_Class_Accuracy(ipretr,inLV(ipretr));
       end
       [~,opt_pretr]=max(Acc_nLV)
       opt_nlv=inLV(opt_pretr)+1

       %plot output

       %rebuild model on entire cal set and apply to test set
       output_PLSDA=Npls_da_N_SEP_pre_SS(Xcal,Ycal,Xtest,Ytest,maxnlv);
       Pred_Ytest=output_PLSDA{1,opt_pretr}.Pred_class_t(:,opt_nlv);


      %plot cross validated accuracies for each pretreatment
        figure, subplot(1,2,1),title('Calibration Set')
        plot(Fold_Mean_Class_Accuracy','-o','LineWidth',2),
        hold on,plot(opt_nlv,Fold_Mean_Class_Accuracy(opt_pretr,opt_nlv),'*k','MarkerSize',20)
        xlabel('Number of Latent Variables','Fontsize',20),
        ylabel('Mean Class Accuracy ','Fontsize',20)
        legend({'none', 'SNV', '1D', '2D', 'SNV + 1D', 'SNV+2D', 'MSC'},'Location','best')
       
        ylim([0,100])
        ax = gca;
        ax.FontSize = 16;
        
        subplot(1,2,2),title('Test Set')
        confusionchart(double(output_PLSDA{1,1}.Yval),Pred_Ytest)


end


%applying PLS-DA with pretreatments: none, SNV, 1D, 2D, SNV + 1D, SNV+2D, MSC
%Applying random cross validation niters times
function [Fold_CV, Fold_SeSp, Fold_Mean_Class_Accuracy,Mean_Class_Accuracy,B_all]=Npls_da_N_SEP_pre_randCV_SS(Xcal,Ycal,maxnlv,niters)
Xcal=double(Xcal);

for i=1:niters
% Set the random seed to 1 for reproducibility
rng(i);
Ind_RP=randperm(size(Xcal,1));
Per30=floor(size(Xcal,1)*0.3);

Xcal2=Xcal(Ind_RP(1:Per30),:);
Ycal2=Ycal(Ind_RP(1:Per30),:);
Xcal1=Xcal(Ind_RP(Per30+1:end),:);
Ycal1=Ycal(Ind_RP(Per30+1:end),:);

output1 = Npls_da_N_SEP_pre_SS(Xcal1,Ycal1,Xcal2,Ycal2,maxnlv);

    for j=1:7
    CC_CV(j,:,i)=output1{1,j}.CC_test;
    SeSp_CV(j,:,i)=output1{1,j}.SeSp;
   for k=1:20
    Mean_Class_Accuracy(j,k,i)=output1{1,j}.SS{1,k}.Perclass_mean_accuracy;
    
   end
   B_all{i}{j}=output1{1,j}.B;
    end
 i   
end

Fold_CV=mean(CC_CV,3);
Fold_SeSp=mean(SeSp_CV,3);
Fold_Mean_Class_Accuracy=mean(Mean_Class_Accuracy,3);

end

%applying PLS-DA with pretreatments: none, SNV, 1D, 2D, SNV + 1D, SNV+2D, MSC
function [output]=Npls_da_N_SEP_pre_SS(Xcal,Ycal,Xval,Yval,maxnlv);
Xcal=double(Xcal);Xval=double(Xval);

[output{1}] = Npls_da_N_SEP_SS(Xcal,Ycal,Xval,Yval,maxnlv);
[output{2}] = Npls_da_N_SEP_SS(SNV(Xcal),Ycal,SNV(Xval),Yval,maxnlv);
[output{3}] = Npls_da_N_SEP_SS(savgol(Xcal,15,3,1),Ycal,savgol(Xval,15,3,1),Yval,maxnlv);
[output{4}] = Npls_da_N_SEP_SS(savgol(Xcal,15,3,2),Ycal,savgol(Xval,15,3,2),Yval,maxnlv);
[output{5}] = Npls_da_N_SEP_SS(savgol(SNV(Xcal),15,3,1),Ycal,savgol(SNV(Xval),15,3,1),Yval,maxnlv);
[output{6}] = Npls_da_N_SEP_SS(savgol(SNV(Xcal),15,3,2),Ycal,savgol(SNV(Xval),15,3,2),Yval,maxnlv);
[output{7}] = Npls_da_N_SEP_SS(MSC(Xcal,mean(Xcal)'),Ycal,MSC(Xval,mean(Xcal)'),Yval,maxnlv);

close all 
end
%% functions

%function to apply PLS-DA on with multiple classes by Aoife Gowen 2016;
%modified by Ana Herrero-Langreo, 2019
%based on NPLS (modelling many Ys at once)
%modified by AG 2020 to include specificity and selectivity
%modified by AG 2026 to include precision, recall, F1

function [output] = Npls_da_N_SEP_SS(Xcal,Ycal,Xval,Yval,maxnlv)
% INPUT
%   Xcal,Ycal,Xval,Yval: spectra (X) and references (Y) for calibration (cal) and
%   validation (val)
%   maxnlv: maximum number of latent variables
%
% [output] = Npls_da_N_SEP(Xcal,Ycal,Xval,Yval,maxnlv)


Ycal_unique=unique(Ycal);


%Labels will be intially replaced by numbers starting from 1.
%Then the Y matrix will be constructed with 0s and 1s for each class
%They will be reconverted into the actual labels at the end

nclasses=size(Ycal_unique,1);

Ycal_num=zeros(size(Ycal,1),1);
Yval_num=zeros(size(Yval,1),1);


for k=1:nclasses
Ycal_cl(:,k)=(Ycal==Ycal_unique(k));
Yval_cl(:,k)=(Yval==Ycal_unique(k));

Ycal_num(Ycal==Ycal_unique(k))=k;
Yval_num(Yval==Ycal_unique(k))=k;
end

B = NplsN(Xcal,Ycal_cl,20);

for j=1:20

Y_pred_cal=Xcal*squeeze(B(:,:,j)) ;
Y_pred_test=Xval*squeeze(B(:,:,j)) ;

%working out thresholds and probabilities here

for k=1:nclasses

Ind1_c=find(Ycal_cl(:,k)==1);
Ind0_c=find(Ycal_cl(:,k)==0);

Y_pred_class1=Y_pred_cal(Ind1_c,k);
Y_pred_class0=Y_pred_cal(Ind0_c,k);

c0_m(j,k)=mean(Y_pred_class0);
c0_sd(j,k)=std(Y_pred_class0);

c1_m(j,k)=mean(Y_pred_class1);
c1_sd(j,k)=std(Y_pred_class1);

%basing on distances
Py_0=1./(sqrt(2*pi)*c0_sd(j,k)) * exp(-0.5*((Y_pred_cal(:,k)-c0_m(j,k))/c0_sd(j,k)).^2); 
Py_1=1./(sqrt(2*pi)*c1_sd(j,k)) * exp(-0.5*((Y_pred_cal(:,k)-c1_m(j,k))/c1_sd(j,k)).^2); 

%proportion of 0's & 1's
P0(k)=size(Ind0_c,1)/(size(Ind1_c,1)+size(Ind0_c,1));
P1(k)=size(Ind1_c,1)/(size(Ind1_c,1)+size(Ind0_c,1));

P0_y=Py_0*P0(k)./(Py_0*P0(k)+Py_1*P1(k));
P1_y=Py_1*P1(k)./(Py_0*P0(k)+Py_1*P1(k));


PO_Y(:,j,k)=P0_y;
P1_Y(:,j,k)=P1_y;
Y_pred_cal_all(:,j,k)=Y_pred_cal(:,k);
Y_pred_test_all(:,j,k)=Y_pred_test(:,k);

%calculate threshold based on quasi normal distribution 

threshold(j,k)=(c1_m(j,k)-2*c1_sd(j,k) - (c0_m(j,k)+2*c0_sd(j,k)))/2+(c0_m(j,k)+2*c0_sd(j,k));

Err_cal(:,k)=((Y_pred_cal(:,k)>threshold(j,k))~=Ycal_cl(:,k));
Err_test(:,k)=((Y_pred_test(:,k)>threshold(j,k))~=Yval_cl(:,k));


end
Err_per_cal(j)=100*sum(Err_cal(:))/numel(Err_cal);
Err_per_test(j)=100*sum(Err_test(:))/numel(Err_test);
end

%figure,plot(Err_per_cal,'o'),xlabel('Number of LV'),ylabel('% Error')

%maxnlv = input('Maximum number of LVs?: ', 's');
%maxnlv=str2double(maxnlv);
%maxnlv=20
for nlv=1:maxnlv
for k=1:nclasses

%taking optimal #LV, predict probablity of each class for each sample
%which class is the most likely for each sample
Py_0_c=1./(sqrt(2*pi)*c0_sd(nlv,k)) * exp(-0.5*((Y_pred_cal_all(:,nlv,k)-c0_m(j,k))/c0_sd(nlv,k)).^2); 
Py_1_c=1./(sqrt(2*pi)*c1_sd(nlv,k)) * exp(-0.5*((Y_pred_cal_all(:,nlv,k)-c1_m(j,k))/c1_sd(nlv,k)).^2); 

P1_y_c(:,k,nlv)=Py_1_c*P1(k)./(Py_0_c*P0(k)+Py_1_c*P1(k));

Py_0_t=1./(sqrt(2*pi)*c0_sd(nlv,k)) * exp(-0.5*((Y_pred_test_all(:,nlv,k)-c0_m(nlv,k))/c0_sd(nlv,k)).^2); 
Py_1_t=1./(sqrt(2*pi)*c1_sd(nlv,k)) * exp(-0.5*((Y_pred_test_all(:,nlv,k)-c1_m(nlv,k))/c1_sd(nlv,k)).^2); 

P1_y_t(:,k,nlv)=Py_1_t*P1(k)./(Py_0_t*P0(k)+Py_1_t*P1(k));

end
end

[Max_prob_c,Pred_class_c]=max(P1_y_c,[],2);
[Max_prob_t,Pred_class_t]=max(P1_y_t,[],2);

cl_error_cal=squeeze(sum(Pred_class_c~=Ycal_num)/size(Ycal,1));
cl_error_test=squeeze(sum(Pred_class_t~=Yval_num)/size(Yval,1));



CC_all_c= 100.*(1-cl_error_cal);
CC_all_t=100.*(1-cl_error_test);

%convert predicted classes to actual class labels
Pred_class_c_actual=zeros(size(Pred_class_c,1),20);
Pred_class_t_actual=zeros(size(Pred_class_t,1),20);
for k=1:nclasses
    for j=1:20
Pred_class_c_actual((Pred_class_c(:,:,j)==k),j)=Ycal_unique(k);
Pred_class_t_actual(Pred_class_t(:,:,j)==k,j)=Ycal_unique(k);
    end
end

    %calculate sensitivity, specificity
if numel(unique(Yval))>1
for j=1:20
[SS{j}]=SensSpecConf_Nclass(Yval,Pred_class_t_actual(:,j));
SeSp(j)=SS{1,j}.Sensitivity(1)*SS{1,j}.Specificity(1);
end
output.SeSp=SeSp;
end

if numel(unique(Yval))>2
for j=1:20
[SS{j}]=SensSpecConf_Nclass(Yval,Pred_class_t_actual(:,j));
for k=1:numel(unique(Yval))
    SeSp(j,k)=SS{1,j}.Sensitivity(k)*SS{1,j}.Specificity(k);
end
output.SeSp(j)=mean(SeSp(j,:));
end
end
    
if numel(unique(Yval))<1
    for j=1:20
    SS{j}=[Yval';Pred_class_t_actual(:,j)'];
    end
end

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
output.nlv=nlv;
output.CC_test=CC_all_t;
output.CC_cal=CC_all_c;
output.Pred_class_c=Pred_class_c_actual;
output.Pred_class_t=Pred_class_t_actual;
output.Ycal=Ycal;
output.Yval=Yval;

end

function [B,T,TT,P,W,Q] = NplsN(X,Y,A)
%NIPALS PLS
v=[];W=[];T=[];TT=[];P=[];Q=[];B=[];
for a = 1:A,
v = X'*Y(:,1);
W(:,a) = v/sqrt(v'*v);
T(:,a) = X*W(:,a);
TT = T(:,a)'*T(:,a);
P(:,a) = X'*T(:,a)/TT;
X = X-T(:,a)*P(:,a)';
Q(a,:) = T(:,a)'*Y/TT;
B(:,:,a) = W*inv(P'*W)*Q;
end
end

function[output]=SensSpecConf_Nclass(Actual,Predicted)
%calculates confusion matrix, sensitivity & specficty from output of PLS-DA model built with Npls_da_N 
%AGowen Feb 2020
%updated AGOWEN 2026 to include precision, F1

nclasses=unique(Actual);


for i=1:size(nclasses,1)
    for j=1:size(nclasses,1)
        Find_equal=(Predicted==nclasses(i));
Confusion_matrix(j,i)=sum(Find_equal(Actual==nclasses(j)));
    end
end

CC=100*sum(diag(Confusion_matrix))/sum(Confusion_matrix(:));

%per class accuracy
for i=1:size(nclasses,1)
    Actual_ind=find(Actual==nclasses(i));
    
Per_class_accuracy(i)=100*sum(Predicted(Actual_ind)==nclasses(i))/numel(Actual_ind);
end

%Calculate Sensitivity & Specificity
for i=1:size(nclasses,1)
     ni2=1:size(nclasses,1);
    ni=setdiff(ni2,i);
   
   TP_class(i) = Confusion_matrix(i,i);
   FN=[];FP=[];TN=[];
   for j=1:size(ni,1)
   FN = [FN,Confusion_matrix(i,ni(j))];
   FP = [FP,Confusion_matrix(ni(j),i)];
   TN = [TN;Confusion_matrix(ni,ni(j))];
   end
   FN_class(i)=sum(FN);
   FP_class(i)=sum(FP);
   TN_class(i)=sum(TN(:));
   
Sensitivity(i) = TP_class(i)/(TP_class(i)+FN_class(i));
%probability test is positive if disease is present 

Specificity(i) = TN_class(i)/(FP_class(i)+TN_class(i));
%probability test is negative if disease not present

Precision(i) = TP_class(i)/(TP_class(i)+FP_class(i));

F1(i) = 2*TP_class(i)/(2*TP_class(i)+FP_class(i)+FN_class(i));

end

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

Xsnv=(X-mean(X,2))./(std(X,0,2));
end

function [z]=savgol(x,F,K,Dn)
% function to calculate n-th order derivative using the Savitzky-Golay coefficients
% [X]=saisir_derivative(X1,window_size,polynom_order,derivative_order)
%
% Input:
%
% x:spectra matrix (n x p)
% F: window_size:(integer) number of data points involved in the calculation (odd number) 
% k: polynom_order:(integer) order of the fitting polynom 
% Dn: derivative_order: (integer, normally 1 or 2) order of the derivative
%
% Output:
% z : transformed data matrix (n rows)
%
% The function assumes that X is a matrix of digitized signals (such as
% spectra) with constant intervals of digitization.
%
% Example:
% z=savgol(x,15,3,2);
% Compute the second derivative using a polynom of power 3 as model
% and a window size of 15
%
% AUTHOR: ANA HERRERO-LANGREO, UCD,2018. Adapted from ALAIN COLLET, 2005 (function from SAISIR environment, written as saisir]=saisir_derivative(saisir1,K,F,Dn))
% modified to be compatible with functions using the eigenvector function savgol.

%
%testing F is odd
aux=F;
gap=floor(F/2);
if(gap==F/2);
    F=F+1;
    xdisp('Warning : the parameter "windowsize" must be odd. It has been set at ',F, ' instead of its initial value ', aux); 
end

[b,g]=sgolaycoef(K,F);
[nrow,ncol]=size(x);
M=zeros(nrow,ncol);
z=zeros(1,ncol);
for i=1:nrow
    
    y=x(i,:);
    for j =(F+1)/2:ncol-(F-1)/2%Calculate the n-th derivative of the i-th spectrum
        if Dn==0
            z(j)=g(:,1)'*y(j - (F+1)/2 + 1: j + (F+1)/2 - 1)';
        else
            z(j)=Dn*g(:,Dn+1)'*y(j - (F+1)/2 + 1:j + (F+1)/2 - 1)';
        end
    end
    M(i,:)=z;
end

z=M;
%saisir=selectcol(saisir,((F+1)/2):(ncol-(F-1)/2));

end

%%%%%%%%%%%%%
function [B,G] = sgolaycoef(k,F)
%sgolaycoef         - Computes the Savitsky-Golay coefficients
%function [B,G] = sgolaycoef(k,F) 
%where the polynomial order is K and the frame size is F (an odd number)
%No direct use

W = eye(F);
s = fliplr(vander(-(F-1)./2:(F-1)./2));
S = s(:,1:k+1);   % Compute the Vandermonde matrix

[Q,R] = qr(sqrt(W)*S,0);

G = S*inv(R)*inv(R)'; % Find the matrix of differentiators

B = G*S'*W; % Compute the projection matrix B

end

function[Xcor]=MSC(X,Target_mean)
Xcal=[ones(size(Target_mean)),Target_mean];
Xcor=zeros(size(X));
for i=1:size(X,1)
Ycal=X(i,:)';   
b = regress(Ycal,Xcal);
Xcor(i,:)=((Ycal-b(1)))/b(2);
end
end