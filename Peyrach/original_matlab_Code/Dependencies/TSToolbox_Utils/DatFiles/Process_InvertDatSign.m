function [returnVar,msg] = Process_InvertDatSign(fbasename,nbChan,varargin)

% USAGE:
%     Process_InvertDatSign(fname)
%     This function inverts the sign of the recording (e.g. for Neuralynx)
%
% INPUTS:
%     fname:        dat file name (with or without '.dat' extension)


% Copyright (C) 2019 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%Parameters
chunk = 3e6; %chunk size, 1e7 needs a lot of RAM & GPU memory (depending on the number of reference channel for the median)

datFile = [fbasename '.dat'];
if ~exist(datFile,'file')
    error([datFile ' does not exist'])
end

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
        dat.Data(idx:(ix+1)*chunk*nbChan) = -m;
    else
        m = dat.Data(idx:end);
        dat.Data(idx:end) = -m;
    end
    
    clear m

end
clear dat
fprintf('\n\n')
