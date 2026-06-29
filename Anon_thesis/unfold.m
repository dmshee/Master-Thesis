function[Xuf]=unfold(X)

[x,y,z]=size(X);
Xuf=reshape(X,x*y,z);