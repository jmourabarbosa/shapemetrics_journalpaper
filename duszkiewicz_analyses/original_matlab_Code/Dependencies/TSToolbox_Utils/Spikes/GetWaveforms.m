function wf = GetWaveforms(filename)


wf = Make_MeanWaveF_FromDat(filename);

wf2 = wf;

for i = 1:length(wf)
    
    wf{i} = double(wf{i});
    
end

SaveAnalysis(pwd,'MeanWaveforms',{wf},{'Waveforms'})


Make_WaveFormFeatures(wf)
