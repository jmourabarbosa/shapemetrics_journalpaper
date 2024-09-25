function [returnVar,msg] = Process_RemoveMedianHighFq(fbasename,nbChan,varargin)

% USAGE:
%     Process_RemoveMedianHighFq(fname,nbChan,refchan,options)
%     This function creates a new dat file including only data from one electrode group (as defined in Neuroscope).
%     It can also substract the median of high pass filtered signals from certain channels
%     (e.g. from a given probe).
%
% INPUTS:
%     fname:        dat file name (with or without '.dat' extension)
%     nbChan:       total number of channels
%     refChan:      vector of electrode indices to compute median

%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'highFc'              high-pass frequency divided by Nyquist frequency 
%                           (2 value vector for bandpass, default low-pass
%                           cut-off: 0.9)
%     'isGPU'               true for GPU computing (for filtering)
%    =========================================================================


% Copyright (C) 2013-18 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%Parameters
chunk = 1e6; %chunk size, 1e7 needs a lot of RAM & GPU memory (depending on the number of reference channel for the median)
refChan = [1:nbChan];
highFc  = [300 8000]/10000; %default, assuming sampling frequency of 20kHz
isGPU = 1; %very significant speed increase with GPU

datFile = [fbasename '.dat'];
if ~exist(datFile,'file')
    error([datFile ' does not exist'])
end

% Parse options
for i = 1:2:length(varargin)
  if ~isa(varargin{i},'char')
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help ' fctName ' ''for details).']);
  end
  switch(lower(varargin{i}))
   
    case 'highfc'
      highFc = varargin{i+1};
      if ~isa(highFc,'numeric') && numel(highFc)>2 && any(highFc>1)
        error(['Incorrect value for property ''highFc'' (type ''help ' fctName ' ''for details).']);
      end  
    case 'isgpu'
      isGPU = varargin{i+1};
      if ~ (isa(isGPU,'numeric') || isa(isGPU,'logical'))
        error(['Incorrect value for property ''isGPU'' (type ''help ' fctName ' ''for details).']);
      end  
    case 'refchan'
      refChan = varargin{i+1};
      if ~ isa(refChan,'numeric')
        error(['Incorrect value for property ''refChan'' (type ''help ' fctName ' ''for details).']);
      end 
  end
end


%% Computing filter options if median is substracted
if length(highFc) == 1
    highFc = [highFc 0.9];
end

[bFilt, aFilt] = butter(3, highFc, 'bandpass');
    
       
infoFile = dir(datFile);

nbChunks = floor(infoFile.bytes/(nbChan*chunk*2));
warning off
if nbChunks==0
    chunk = infoFile.bytes/(nbChan*2);
end

dat = memmapfile([fbasename '.dat'],'Format','int16','writable',true);

for ix=0:nbChunks
   fprintf('.')
    % load data in a memory map
    idx = ix*chunk*nbChan+1;
    if ix<nbChunks
        m = dat.Data(idx:(ix+1)*chunk*nbChan);
    else
        m = dat.Data(idx:end);
        chunk = infoFile.bytes/(2*nbChan)-nbChunks*chunk;
        datF = reshape(m,[nbChan chunk]);
    end
    
    datF = reshape(m,[nbChan chunk]);
    datF = double(datF);
    
    if isGPU
        datF = gpuArray(datF);
    end
    
    datF = filter(bFilt, aFilt, datF,[],2);
    datF = fliplr(datF);
    datF = filter(bFilt, aFilt, datF,[],2);
    datF = fliplr(datF);
    refF = median(datF(refChan,:));
    datF = bsxfun(@minus, datF, refF);
    
    if isGPU
        datF = gather(datF);
    end
    
    if ix<nbChunks
        dat.Data(idx:(ix+1)*chunk*nbChan) = datF(:);
    else
        dat.Data(idx:end) = datF(:);
    end
    % write data to disk
    clear m tmpDat refDat datF

end
clear dat
fprintf('\n\n')
