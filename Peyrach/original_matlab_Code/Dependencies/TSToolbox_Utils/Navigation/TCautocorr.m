function ac = TCautocorr(h)

if length(h) > 360 % assuming 360 bins
    h = h(1:360);
end

nBins = length(h);



ac = nan(nBins,1);

for nR = 0:nBins-1 
    hs = circshift(h,nR);
    ac(nR+1) = corr(h,hs);
end