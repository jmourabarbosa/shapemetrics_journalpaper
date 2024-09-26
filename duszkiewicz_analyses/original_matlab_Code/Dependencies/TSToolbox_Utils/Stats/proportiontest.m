function [z,p] = proportiontest(p,p0,n)

z = (p-p0) / sqrt(p0*(1-p0)/n);
[~,p] = ztest(z,0,1);
