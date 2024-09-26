function meanWaveF = Make_MeanWaveF(fbasename,varargin)

% meanWaveF = Make_MeanWaveF(fbasename,option)
% Computes average waveform from spk files
%
% options:
%    'spk' (default) | 'dat': reads spike from the spk file of from the dat file
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        %Consder moving this line to the loop if memory issues (running on
        %single clu index)
        wav = LoadSpikeWaveF([fbasename '.spk.' num2str(ch)], ne,ns,find(clu>1));
        clu = clu(clu>1);

        for c=1:length(nClu)
            %wav = LoadSpikeWaveF([fbasename '.spk.' num2str(ch)], ne,ns,find(clu==nClu(c)));
            w = [];
            for e=1:ne
                %w = [w;mean(squeeze(wav(e,:,:)),2)'];
                w = [w;mean(squeeze(wav(e,:,clu==nClu(c))),2)'];
            end

            %clear wav

            meanWaveF{end+1} = w;
            
            if c>1 && c<11
                fprintf(repmat('\b',[1 length(strg)-1]));
            elseif c>10
                fprintf(repmat('\b',[1 length(strg)]));                
            end
            strg = 'Cell %i done';
            fprintf(strg,c)
        end

    else
        warning([fname ' doesn''t exist']);
    end

    clear spk;

end