function ac = TCcrosscorr(h1,h2)

if rem(h1,2) == 1 % if odd number of bins then make even
    h1 = h1(1:end-1);
end

if rem(h2,2) == 1 
    h2 = h2(1:end-1);
end

nBins = length(h1);
ac = nan(nBins,1);

for nR = 0:nBins-1
    hs = circshift(h2,nR);
    ac(nR+1) = corr(h1,hs,'Type','Pearson');
end