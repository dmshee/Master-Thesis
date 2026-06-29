function [Xcor] = MSC_reworked(X, Target_mean)
% X: nSamples x nBands
% Target_mean: 1 x nBands or nBands x 1

X = double(X);
Target_mean = double(Target_mean(:));   % nBands x 1

[nSamples, nBands] = size(X);
if numel(Target_mean) ~= nBands
    error('MSC_reworked:SizeMismatch', ...
        'Target_mean length (%d) must equal number of bands in X (%d).', ...
        numel(Target_mean), nBands);
end

Xcal = [ones(nBands,1), Target_mean];   % nBands x 2
Xcor = zeros(nSamples, nBands);

for i = 1:nSamples
    Ycal = X(i,:).';                     % nBands x 1
    b = Xcal \ Ycal;                     % least-squares fit
    Xcor(i,:) = ((Ycal - b(1)) / b(2)).';
end
end