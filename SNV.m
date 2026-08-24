%% SNV function to carry out Standard Normal Variate preprocessing on a group of spectra, X
% X should be a 2D matrix with one spectrum per row.
% written by A. Gowen

function [Xcor] = SNV(X)
X = double(X);
[r, c] = size(X);
X_rowmean = mean(X, 2);
X_rowstd = std(X, 0, 2);
X_rowstd(X_rowstd == 0) = 1;

Xcor = (X - repmat(X_rowmean, 1, c)) ./ repmat(X_rowstd, 1, c);
end

