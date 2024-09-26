%  USAGE
%
%    data = Make_MeanWaveF_FromDat(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nSamples'    number of samples per waveform (default: 40)
%     'peakIndex'   index of the peak (default: 16)
%    =========================================================================

% Copyright (C) 2019 by Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

  
function meanWaveF = Make_MeanWaveF_FromDat(fbasename,varargin)

%Default values:
nsamples = 40;
peakindex = 16;

if nargin < 1 | mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help Make_MeanWaeveF_FromDat'' for details).');
end

% Parse options
for i = 1:2:length(varargin)
  if ~isa(varargin{i},'char')
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help Make_MeanWaeveF_FromDat'' for details).']);
  end
  switch(lower(varargin{i}))
    case 'nsamples'
      nsamples = varargin{i+1};
      if ~isa(duration,'numeric') | length(nsamples) ~= 1 | nsamples < 1
        error('Incorrect value for property ''nsamples'' (type ''help Make_MeanWaeveF_FromDat'' for details).');
      end
    case 'peakindex'
      peakindex = varargin{i+1};
      if ~isa(frequency,'numeric') | length(peakindex) ~= 1 | peakindex > nsamples 
        error('Incorrect value for property ''peakindex'' (type ''help Make_MeanWaeveF_FromDat'' for details).');
      end
     otherwise
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinary'' for details).']);
  end
end


xml_data = LoadXml(fbasename);

shankIx = [1:length(xml_data.SpkGrps)];
meanWaveF = {};

for ch=shankIx

    ne = length(xml_data.SpkGrps(ch).Channels);
    ns = xml_data.SpkGrps(ch).nSamples;

    disp(['Shank #' num2str(ch)])

    fname = [fbasename '.clu.' num2str(ch)];

    if exist(fname,'file')

        clu = load(fname);
        clu = clu(2:end);
        nClu = unique(clu(clu>1));
        
        clu = clu(clu>1);

        for c=1:length(nClu)
            wav = load_spk_from_dat_wavef(fbasename, ch,nClu(c),nsamples,peakindex,1);
            wav = squeeze(mean(wav,3));
            meanWaveF{end+1} = wav;
        end

    else
        warning([fname ' doesn''t exist']);
    end

end