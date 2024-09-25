function MakeSpikeWidth(varargin)

% MakeSpikeWidth - computes and saves cell waveform features : half peak width and peak to peak
% both of them in ms.
%
%  USAGE
%
%    MakeSpikeWidth(fbasename)
%
%    fbasename (optional)   file base name (default: folder's name)
%

% Copyright (C) 20012-2018 by Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


if isempty(varargin)
    [~,fbasename,~] = fileparts(pwd);
else
    fbasename = varargin{1};
end

[spkWidth,pk2pk,halfPk,meanWaveF,maxIx] = Make_MeanWaveF(fbasename);
         
info = {'Peak 2 Peak';'spike Width (inverse of peak spectrum)';'half peak';'mean wave forms';'maximal spk channel'};
SaveAnalysis(pwd,'SpikeWaveF',{pk2pk; spkWidth;halfPk;meanWaveF;maxIx},{'pk2pk'; 'spkWidth'; 'halfPk'; 'meanWaveF';'maxIx'},info);
