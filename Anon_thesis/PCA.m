function[L,S,Ev]=PCA(X)
%function to carry out PCA on a matrix of spectra (X) by A.Gowen
%Inputs: X = Spectra (or unfolded spectral image)
%Outputs: L = PCA loadings, S = PCA scores, Ev = PCA eigenvalues

X=double(X);
[r,c]=size(X);

%mean centre X
X_mean=mean(X);
X_mncn=X-repmat(mean(X),r,1);

%calculate covariance matrix
Xcov=cov(X_mncn);

%calculate L & Ev by SVD
[U,Sev,V]=svd(Xcov);
L=V;
%calculate S
S=X_mncn*L;
Ev=diag(Sev);