function MakeSpikeData_FromKlusters(varargin)

% MakeSpikeData_FromKlusters(filebasname)
% construct a tsdArray of spike timing from the Neurosuite feils and save it in Analysis/SpikeData.mat
% 
% INPUTS:
%     filebasename (optional): file base name (folder name is default)
                               
% Adrien Peyrache, 2011

if isempty(varargin)
    [~,fbasename,~] = fileparts(pwd);
else
    fbasename = varargin{1};
end

% Default values
shanks = [];
destFolder = pwd;

if nargin>1 & mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters.');
end
% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property.']);
  end
  switch(lower(varargin{i})),
    case 'shanks',
      shanks = varargin{i+1};
      if ~isa(shanks,'numeric')
        error('Incorrect value for property ''shanks''.');
      end
    case 'destfolder'
      destFolder = varargin{i+1};
      if ~isa(destFolder,'char')
        error('Incorrect value for property ''destFolder''.');
      end
      if ~strcmp(destFolder(end),'/')
        destFolder = [destFolder '/'];
        destFolder = [destFolder fbasename];   
      end
  end
end

if ~exist(destFolder,'dir')
    mkdir(destFolder)
end
if isempty(shanks)
    info = LoadXml([fbasename '.xml']);
    nbSh = length(info.ElecGp);
    shanks = (1:nbSh);
end

%muaFile = [analysisDir filesep 'MUAData.mat'];

[S,shank,cellIx ] = LoadSpikeData2(fbasename,shanks);

fprintf('Number of Cells: %d\n',length(S))
SaveAnalysis(destFolder,'SpikeData',{S,shank,cellIx},{'S','shank','cellIx'});
%save(muaFile,'MUA');