function z = savgol(x, F, K, Dn)
%SAVGOL Compute Savitzky-Golay smoothing/derivative for row-wise spectra.
%   z = savgol(x, F, K, Dn)
%   x  : spectra matrix (nSpectra x nVariables)
%   F  : odd window size
%   K  : polynomial order
%   Dn : derivative order (0, 1, 2, ...)

if nargin < 4 || isempty(Dn)
    Dn = 0;
end

x = double(x);

if mod(F, 2) == 0
    originalF = F;
    F = F + 1;
    warning('savgol:EvenWindowAdjusted', ...
        'Window size must be odd. Using F = %d instead of %d.', F, originalF);
end

[~, g] = sgolaycoef(K, F);
[nrow, ncol] = size(x);
M = zeros(nrow, ncol);

halfWindow = (F - 1) / 2;
leftPad = repmat(x(:, 1), 1, halfWindow);
rightPad = repmat(x(:, end), 1, halfWindow);
xPad = [leftPad, x, rightPad];

for i = 1:nrow
    y = xPad(i, :);

    for j = 1:ncol
        idx = j : j + F - 1;
        if Dn == 0
            M(i, j) = g(:, 1)' * y(idx)';
        else
            M(i, j) = Dn * g(:, Dn + 1)' * y(idx)';
        end
    end
end

z = M;
end

function [B, G] = sgolaycoef(k, F)
%SGOLAYCOEF Compute Savitzky-Golay projection and differentiator matrices.

W = eye(F);
s = fliplr(vander(-(F - 1) / 2 : (F - 1) / 2));
S = s(:, 1:k + 1);

[~, R] = qr(sqrt(W) * S, 0);

G = (R \ (R' \ S'))';
B = G * S' * W;
end