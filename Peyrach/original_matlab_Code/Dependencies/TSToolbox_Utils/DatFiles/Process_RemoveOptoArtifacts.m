function Process_RemoveOptoArtifacts(datName,nchannels,optoTimes,maxT)

% Process_RemoveOptoArtifacts
% 
%  USAGE
%
%    Process_RemoveOptoArtifacts(datName,nchannels,optoTimes)
%
%    datName        file base name of the dat file (e.g. 'datName' for datName.dat).
%    nchannels      number of channels
%    optoTimes      2 column matrix of start and end times of opto stims (in number of samples)
%    maxT           duration of artifact


% Copyright (C) 2018 by Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

optoDur = diff(optoTimes,1,2);
N       = size(optoTimes,1);
borderFilter = gausswin(101,5);
borderFilter = [borderFilter(1:50);ones(maxT,1);borderFilter(51:end)];

dat = memmapfile(datName,'Format','int16','writable',true);
    
for ii = 1:nchannels
    
   datCh    = dat.Data(ii:nchannels:end);
   
   onIx         = optoTimes(:,1) + repmat((-100:maxT),[N 1]);
   onTemplate   = mean(datCh(onIx));
   onTemplate   = onTemplate(:).*borderFilter;
   onIx         = onIx';
   datCh(onIx(:))  = datCh(onIx(:)) - int16(repmat(onTemplate,[N 1]));
   
   onIx         = optoTimes(:,2) + repmat((-100:maxT),[N 1]);
   onTemplate   = mean(datCh(onIx));
   onTemplate   = onTemplate(:).*borderFilter;
   onIx         = onIx';
   datCh(onIx(:))  = datCh(onIx(:)) - int16(repmat(onTemplate,[N 1]));
   
   dat.Data(ii:nchannels:end) = datCh;
end