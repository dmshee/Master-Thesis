%% SNV function to carry out Standard Normal Variate preatment on a group of spectra, X
% X should be an 2D matrix
% written by A. Gowen

function[Xcor]=SNV(X)
[r,c]=size(X)
X_rowmean=mean(X,2);
X_rowstd=std(X,[],2);

Xcor=(X-repmat(X_rowmean,1,c))./repmat(X_rowstd,1,c);

