function SpatBits=SpatialInfo_1D(pf,occ)

%  USAGE
%  	SpatBits=SpatialInfo(pf,occ)
%  	
%  compute for each cells the number of spatial bits per spikes
%  
%  INPUT:
%  	pf: placeField (In Hz per pixel), organized oolumnwise if multiple
%  	occ: Occupancy probability
%  OUTPUT:
%  	SpatBits: the spatial information (in bits/spikes)
%  
%  Adrien Peyrache 2021

occ = occ(:);
occ = occ/sum(occ);
f = occ'*pf;
pf = pf./repmat(f,[size(pf,1) 1]);
pf(pf==0) = NaN;
SB = (repmat(occ,[1 size(pf,2)]).*pf).*log2(pf);
SpatBits = nansum(SB);

