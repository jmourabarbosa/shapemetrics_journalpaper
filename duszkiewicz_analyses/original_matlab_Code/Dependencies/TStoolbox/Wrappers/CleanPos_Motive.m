% Dependency of LoadPosition.
% [CleanWhl GoodRanges] = CleanWhlFR(pos, StretchLen, JumpSize,Gap);
%
% "cleans up" a position file by interpolating missing stretches
% up to StretchLen long (default 20), for which the endpoints 
% don't differ by more than JumpSize, (default 30).
%
% also returns the ranges where the whl file is valid (in .whl units)
% GoodRanges which gives start and end samples of the good ranges
% (so pos(GoodRanges) has no NaN values).
%
% if there are any very high derivatives left over, it warns you


% Before useing interp1, remove the data which has (-1,-1) or big Gap between two continous rows.
% Assing the (-Gap, -Gap) to the row which has big Gap or (-1,-1).
% interporate.
% if the Gap is longer than StretchLen in terms of number of .whl rows or longer than JumpSize in terms of distance, remove the interporated values.


% changed by A Peyrache, 2012-2017

function [cWhl, GoodRanges] = CleanPos_Motive(Whl, StretchLen, JumpSize, Gap)

% If the gap between the good strech is more than StrethcLen in terms of Whl row number,remove interporated values.
if nargin<2
	StretchLen = 30;
end

% If the Gap between the good strech is more than JumpSize, remove interporated values.
if nargin<3
    JumpSize = 1;
end

% if the distance between the two contimous rows are more than Gap centimeter, It's a big jump and do not use as an input for inpterp1.
if nargin<4
	Gap = 1;
end

nWhl = size(Whl,1);

% interpolate missing values or large jumps.
% I should use distance, not the one dimentinal projection of trajectory, by the way.

whltemp = Whl;
dist = sqrt(diff(whltemp(:,1)).^2+diff(whltemp(:,2)).^2);
BigJump = dist>Gap;

Good = find(~([BigJump;0] | [0;BigJump]));
Bad = find(([BigJump;0] | [0;BigJump]));

whltemp(Bad,1:2) = -Gap;

% Give -1 outside of the interpolation.

if length(Good)<2 
    cWhl(:,1:2) = NaN(size(Whl,1),2);
else
	cWhl(:,1:2) = interp1(Good, Whl(Good,1:2), 1:nWhl, 'linear', -1);

end

% find missing stretches
dGood = [-(whltemp(1,1)==-Gap) ; diff(whltemp(:,1)>-Gap)];
BadStart = find(dGood<0);
BadEnd = find(dGood>0)-1;
% if last point is bad, need to finish off BadEnd
if Whl(end,1)==-1
	BadEnd = [BadEnd; nWhl];
end

if length(BadStart)>length(BadEnd)
	BadEnd = [BadEnd; nWhl];
end


% find ranges to chuck
% jump size ...
if any(BadStart>0)

    StartInd = clip(BadStart-1, 1, nWhl); % StartInd and EndInd give the 
    EndInd = clip(BadEnd+1, 1, nWhl);     % points you are interpolating between
    
	dist = sqrt((Whl(StartInd,1)-Whl(EndInd,1)).^2+(Whl(StartInd,2)-Whl(EndInd,2)).^2);
    ToChuck = find(BadEnd-BadStart>=StretchLen | dist > JumpSize);
	% chuck em

	for i=ToChuck(:)'
       	cWhl(BadStart(i):BadEnd(i),1:2) = NaN;
	end
end

dcGood = diff([0; ~isnan(cWhl(:,1)); 0]);
GoodStart = find(dcGood>0);
GoodEnd = find(dcGood<0)-1;
GoodRanges = [GoodStart, GoodEnd];


return
