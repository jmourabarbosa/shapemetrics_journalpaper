function [h,B,mu,h0] = HeadDirectionField_BD(tsa,ang,GoodRanges,varargin)

% USAGE
%     [h,d,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges)
%     
%     Inputs:
%     tsa: a ts object (typically a cell!)
%     ang: a tsd object of head orientation
%     GoodRanges: an intervalSet object defining the time of valid pos tracking
%
%     options:
%     [h,d,mu,p,kappa] = HeadDirectionField(tsa,ang,GoodRanges,nbBins,sdSmooth)
%     nbBins: number of angle bins (default=360)
%     sdSmooth: s.d. of tuning curve smoothing (in number of bins, default=6)

% This variant of the function calculates a bidirectional tuning curve
% following Jacob et al 2015

%Default values
nbBins = 360;
sdSmooth = 3; %s.d. of smoothing filter in number of bins

h0 = [];
if ~isempty(varargin)
    if isnumeric(varargin{1}) && length(varargin{1})==1
        nbBins=varargin{1};
    else
        error('Nb of bins must be numeric')
    end
    if length(varargin) == 2
        if isnumeric(varargin{2}) && length(varargin{2})==1
            sdSmooth=varargin{2};
        else	
            error('S.d. of smoothing filter must be numeric')
        end
    end
end

tsa = Restrict(tsa,GoodRanges);
B = 2*pi*(0:1:nbBins-1)'/nbBins;

if ~isempty(Range(tsa))
    
    ang = Restrict(ang,GoodRanges);
    angt = Restrict(ang,tsa);
    angrange = Range(ang);
    ang = mod(Data(ang),2*pi);
    angt = mod(Data(angt),2*pi);

    
    h0 = hist(mod(ang*2,2*pi),B);
    h = hist(mod(angt*2,2*pi),B);
    
    dt = 1/median(diff(angrange));
    h = dt*h./h0;
    
    h(isnan(h))=0;
    h = h(:);

    if sdSmooth
        h = gaussFiltAng(h,sdSmooth,1);
    end
    h = [h;h(1)];
    B = [B;B(1)];

    [mu, kappa, p] = CircularMean(angt);
else
    h = zeros(nbBins+1,1);
    B = [B;B(1)];
    mu = NaN;
    kappa = NaN;
    p = NaN;
end