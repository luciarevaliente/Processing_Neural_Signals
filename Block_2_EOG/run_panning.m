clear all; clc
realtime=0; %flag to select between processing recording in real time, or reading in audio "sample" instead
time=0; %in seconds, duration of realtime processing
algo=db_panning; 
algo.auto = 1;

initdata.BufferSize=512;

[aPR] = init_aPR(initdata.BufferSize); %initialize audioplayerrecorder (M-Audio Drivers need to be installed)

if realtime
    [recdata,playdata] = realtime_processing(aPR,time,algo,initdata,[], [1 2],[1 2]);
else
    sample=audioread('AIP_song.wav');
    [playdata] = realtime_sample_processing(aPR,sample, algo,initdata,[], [1 2]);
end

