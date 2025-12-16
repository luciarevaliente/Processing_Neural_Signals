%% Task 12 – Run Delay-and-Sum Beamformer
clear; clc;

% --- Initialization ---
BufferSize = 512;
initdata.BufferSize = BufferSize;

aPR = init_aPR(BufferSize);

% --- Load microphone array recording (Task 11) ---
load('task11_mic_recording.mat');  % mic1, mic2, fs

micRec = [mic1 mic2];

% --- Create beamformer plugin ---
algo = DSB;
algo.MicSpacing = 0.08;   % meters – set to your actual spacing!

% --- Run offline simulated real-time processing ---
playdata = offline_sample_processing( ...
    aPR, micRec, algo, initdata, [], [1 2]);

% --- Cleanup ---
release(aPR);
