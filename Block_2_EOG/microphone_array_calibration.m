%% Task 10: Microphone Array Calibration

% 1. Init the audio device
aPR = init_aPR(); 

% Settings for recording
fs = 44100;             % sample rate (Task 8: 44.1 kHz)
recDuration = 5;        % Recording Duration
frameSize = 512;        % Standard Buffers Size
numFrames = ceil((recDuration * fs) / frameSize);

disp('Switch Calibrator on and put on the microphone.');
disp('Recording starts in 2 seconds...');
pause(2);

%% 2. Aufnahme des Testsignals [cite: 681]
disp('Recording...');

recSignal = []; % Save the signal

% Puffer for the recording
silenceBuf = zeros(frameSize, 2); 

for i = 1:numFrames
    % aPR returns the recording
    [recBuf] = aPR(silenceBuf);
    
    % Attach the puffer
    recSignal = [recSignal; recBuf];
end

disp('Recording finished.');

% Release aPR
release(aPR);

%% 3.Calculation of RMS and Scaling factor
% Assuming calibration of mics 

% Calculate RMS for each channel
rms_ch1 = rms(recSignal(:, 1));
rms_ch2 = rms(recSignal(:, 2));

fprintf('RMS Channel 1: %.4f\n', rms_ch1);
fprintf('RMS Channel 2: %.4f\n', rms_ch2);

% Check: Is there a signal?
if rms_ch1 < 0.001 || rms_ch2 < 0.001
    warning('RMS very low. Record again');
end

% Calculate scaling factor (mic_scale)
scale_factor_1 = 1 / rms_ch1;
scale_factor_2 = 1 / rms_ch2;

mic_scale = [scale_factor_1; scale_factor_2];

fprintf('Calculate scaling factors:\n  Channel 1: %.4f\n  Channel 2: %.4f\n', ...
    mic_scale(1), mic_scale(2));

%% 4. Saving of Factors
save('mic_scale.mat', 'mic_scale');
disp('Saved scaling factord in mic_scale.mat.');