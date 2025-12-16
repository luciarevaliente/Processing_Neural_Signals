%% Verification of Task 9
%% 1. Load Audio File
% Load one of the provided audio files (Task 9, Step 3) [cite: 591]
% Replace 'speech_sample.wav' with the actual filename provided in your course folder.
filename = 'AIP_Song.wav'; 

if exist(filename, 'file')
    [input_signal, fs] = audioread(filename);
else
    error('File not found. Please ensure the audio file is in your specific path.');
end

%% 2. Initialize the Plugin and Set Parameter
% Create an instance of your plugin
plugin = dB_panning;

% Set a test gain value (e.g., -3 dB) to verify the logic
% According to your code: Negative slider -> Left louder (+), Right softer (-)
test_gain = -3; 
plugin.Gain = test_gain;

%% 3. Process the Signal
% We call the process method directly to verify the algorithm's math (Task 9, Step 2).
% Note: For the listening test later, you would use 'offline_sample_processing'.
output_signal = plugin.process(input_signal);

%% 4. Calculate RMS and Verify (Task 9, Step 2)
% Calculate RMS (Root Mean Square) for input (r_rms) and output (p_rms) [cite: 585]
% We verify each channel separately because they should have opposite gains.
r_rms_left  = rms(input_signal(:, 1));
r_rms_right = rms(input_signal(:, 2));

p_rms_left  = rms(output_signal(:, 1));
p_rms_right = rms(output_signal(:, 2));

% Calculate the measured gain in dB: 20 * log10(p_rms / r_rms) [cite: 585]
measured_dB_left  = 20 * log10(p_rms_left / r_rms_left);
measured_dB_right = 20 * log10(p_rms_right / r_rms_right);

%% 5. Display Results
fprintf('--- Verification Results for Slider Gain: %.2f dB ---\n', test_gain);

% Left Channel: Expected to be INCREASED by slider value (if negative)
% Logic: gainLeft = -(-3) = +3 dB
fprintf('Left Channel:\n');
fprintf('  Expected: %+.2f dB\n', -test_gain);
fprintf('  Measured: %+.2f dB\n', measured_dB_left);

% Right Channel: Expected to be DECREASED by slider value
% Logic: gainRight = -3 dB
fprintf('Right Channel:\n');
fprintf('  Expected: %+.2f dB\n', test_gain);
fprintf('  Measured: %+.2f dB\n', measured_dB_right);

%% 6. Listening Test (Task 9, Step 3 & 4)
% To listen to the result as requested in Step 4[cite: 598]:
% sound(output_signal, fs);
