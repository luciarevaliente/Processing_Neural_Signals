function[output,state] = RT_EOG(cmd,state)

% RT_EOG
% To be adapted by students of the Ringpraktikum Neuro-Signale course. When
% called, this function should retreive the next availabe EOG acquisition
% buffer from the Biopac MP36 system, process that buffer and estimate the
% angle of the users gaze in real-time. The function is designed to be
% called repeatedly within a real-time loop of EOG_RT_Framework.m.
%
% Input/Output Variables:
%
% cmd:      Input string used to choose which phase of the function is excuted.
%           'init' executes the intialization phase of the function used to
%           compute initial conditions of variables, filter coefficients,
%           etc that should be computed before the real-time loop commences.
%           'process' executes the main processing phases that should be
%           called within the real-time processing loop
%
% state:    Input/Output structure used to store any variables that need to
%           be stored and passed between subsequent calls of RT_EOG. You
%           are free to add any fields to this structure as you see fit,
%           which can be initialized and set within this function. However,
%           the following fields must be included and must be initiliazed
%           OUTSIDE this function.
%
%           state.fs: Sampling rate of EOG acquisition with BIOPAC.
%
%           state.Nacq_buff: Length of EOG acquisition buffer in samples.
%
%           state.is_online: If 1, the function retreives buffers from the
%           BIOPAC system. If 0, function loads in a pre-recorded EOG
%           signal and processes it as if it was being read from the BIOPAC
%
%           state.fr_idx: Tracks how many buffers have been retrieved and
%           processed since the start of execution. Can be used to avoid
%           computing gaze angles while the EOG signal is fluctuating due
%           to system switch-on artefacts.
%
% output:   Output structure comprising an variable you would like to have
%           accessible once the function returns. You can add any fields to
%           this structure you desired, however the following must be
%           included for correct operation of with EOG_RT_Framework
%
%           output.V_raw: the raw EOG data retrieved from the acquistion
%           buffer
%
%           output.V_LP: lowpass filtered version of EOG buffer
%
%           output.V_BP: bandpass version of EOG buffer
%
%           output.edge_idx: vector that is same length as EOG buffer that
%           has contains +1 and -1 on samples where a saccade start and end
%           have respectively been detected, and 0 on others.
%
%           output.ang_est: (scaler) the estimated angle of the user's gaze
%           as of the end of the EOG buffer

persistent b_LP a_LP zi_LP    % Main lowpass filter (35 Hz)
persistent b_baseline a_baseline zi_baseline  % Baseline removal (0.16Hz)
    
if strcmpi(cmd,'process')
    %% Read in next available EOG buffer
    [~, output.V_raw] = biopacAPI(state.is_online,'receiveMPData',state.Nacq_buff); % pull eog data frame
    
    %% Task 3: Apply filtering to Acq Buffer
    % Answer:
        % Problem: The bandpass filtered signal may contain drift or low-frequency noise, especially when eye movements are infrequent.
        % Cause: The baseline removal filter (low-pass at 0.16 Hz) may not fully eliminate slow drift from the signal, leading to misalignment 
        % between the gaze angle and the beamformer's orientation.
    [V_LP_temp, zi_LP] = filter(b_LP, a_LP, output.V_raw, zi_LP);

    % Apply baseline removal filter (extracts slow drift)
    [V_baseline, zi_baseline] = filter(b_baseline, a_baseline, V_LP_temp, zi_baseline);

    % Compute bandpass signal by subtraction (highpass effect)
    output.V_LP = V_LP_temp;
    output.V_BP = V_LP_temp - V_baseline;

    % ----- Global sample indices for this buffer -----
    Nbuff        = state.Nacq_buff;
    fr_idx       = state.fr_idx;
    global_start = (fr_idx-1)*Nbuff + 1;
    global_idx_vec = (global_start : global_start + Nbuff - 1)';
    
    % ----- Append low-pass samples to history for Task 5 averaging -----
    state.pot.hist_values = [state.pot.hist_values; output.V_LP(:)];
    state.pot.hist_index  = [state.pot.hist_index;  global_idx_vec];
    
    % Keep only the most recent max_hist samples
    if numel(state.pot.hist_values) > state.pot.max_hist
        excess = numel(state.pot.hist_values) - state.pot.max_hist;
        state.pot.hist_values(1:excess) = [];
        state.pot.hist_index(1:excess)  = [];
    end
 
    %% Task 4: The saccade edge detection algorithm
    % Answer: 
        % Short filter length causes insufficient smoothing, leading to noisy edge detection and false positives due to high-frequency fluctuations or spikes.
        % High threshold (e.g., 0.01) detects only strong saccades, potentially missing smaller saccades.
        % Low threshold (e.g., 0.0001) makes detection too sensitive, leading to false positives from noise or small changes.
    fr_idx = state.fr_idx;

    % ensure column vectors for processing
    V = output.V_LP(:);
    N = length(V);

    % init edge_idx with same size as V_LP
    edge_idx_col = zeros(N,1);

    if fr_idx <= state.saccade.skip_buffers
        edge_idx_col(:) = 0;
        state.saccade.prev_sample = V(end);
    else
        % 1. derivative dV/dn
        dVLP = zeros(N,1);
        dVLP(1)     = V(1) - state.saccade.prev_sample;
        dVLP(2:end) = diff(V);
        state.saccade.prev_sample = V(end);

        % 2. moving average smoothing
        [ma_output, state.saccade.filter_state] = filter( ...
            state.saccade.ma_coeffs, 1, dVLP, state.saccade.filter_state);

        % 3. threshold on |dV|
        binary_seq = abs(ma_output) > state.saccade.deriv_thresh;  % Nx1 logical

        % 4. delay-and-subtract
        delayed_binary = [state.saccade.prev_binary; binary_seq(1:end-1)];
        edge_idx_col   = double(binary_seq) - double(delayed_binary);

        % 5. save for next iteration
        state.saccade.prev_binary = binary_seq(end);
    end

    % reshape edge_idx back to original shape of V_LP
    output.edge_idx = reshape(edge_idx_col, size(output.V_LP));


    %% Task 5: Estimating the potential changes
    % Append new LP filtered values and their global indices to the history buffer.
    % This allows us to "look back" in time to calculate averages (V_initial/V_final).
    state.pot.hist_values = [state.pot.hist_values; output.V_LP(:)];
    state.pot.hist_index  = [state.pot.hist_index;  global_idx_vec];

    % Maintain a sliding window of fixed size (Circular Buffer logic)
    if numel(state.pot.hist_values) > state.pot.max_hist
        excess = numel(state.pot.hist_values) - state.pot.max_hist;
        state.pot.hist_values(1:excess) = [];
        state.pot.hist_index(1:excess)  = [];
    end
    
    % Initialize local estimate with the current state estimate
    current_V_est = state.pot.V_est; 
    
    % Iterate through every sample in the current block
    for i = 1:N
        curr_global_idx = global_idx_vec(i);
        edge_val = output.edge_idx(i);
        
        % -----------------------------------------------------------
        % 1. TIMEOUT CHECK (Normal Saccade Confirmation)
        % -----------------------------------------------------------
        % If we are 'waiting', a saccade has recently ended (edge -1), and we 
        % are waiting for the eye to settle (handling overshoot/glissade) 
        % to measure the new stable voltage.
        if state.pot.is_waiting
            state.pot.timer = state.pot.timer + 1;

            % Check if enough time has passed (Sloppy period + Averaging window)
            if state.pot.timer >= (state.pot.N_sloppy + state.pot.N_avg)
                % -- Time to measure V_final (The new Peak) --
                % Define the window: Start AFTER the sloppy period to avoid artifacts
                idx_start = state.pot.pending_end_idx + state.pot.N_sloppy;
                idx_end   = idx_start + state.pot.N_avg;
                
                % Extract values from history based on calculated indices
                mask = (state.pot.hist_index >= idx_start) & (state.pot.hist_index <= idx_end);
                vals = state.pot.hist_values(mask);
                
                if ~isempty(vals)
                    V_final = mean(vals);
                    
                    % Smart Logic that has drift correction
                    if isnan(state.pot.last_peak_level)
                        % Case: Very first saccade. 
                        % We have no previous history, so we rely on the immediate Initial vs Final.
                        diff_val = V_final - state.pot.pending_V_initial;
                    else
                        % Case: Subsequent saccades.
                        % instead of comparing (V_final - V_current_initial), we compare
                        % (V_final - V_previous_final).
                        % WHY? Because the signal might have drifted slowly while the eye
                        % was still. By comparing peak-to-peak, we ignore that slow drift.
                        diff_val = V_final - state.pot.last_peak_level;
                    end
                    
                    % Update the total estimated potential change
                    current_V_est = current_V_est + diff_val;
                    
                    % Save this V_final as the stable reference for the 
                    % NEXT movement
                    state.pot.last_peak_level = V_final;
                end
                
                % Reset waiting state (Saccade processing complete)
                state.pot.is_waiting = false;
                state.pot.timer = 0;
            end
        end
        
        % -----------------------------------------------------------
        % 2. SACCADE START (+1)
        % -----------------------------------------------------------
        if edge_val == 1
            % A new saccade is attempting to start.
        
            % Rapid Saccade
            % If 'is_waiting' is true, it means a NEW saccade started before the 
            % old one finished settling.
            if state.pot.is_waiting
                if state.pot.timer < state.pot.N_sloppy
                    % Case A: Too fast (e.g., < 30ms). Likely noise or a blink artifact.
                    % Abort the pending measurement.
                    state.pot.is_waiting = false;
                    
                elseif state.pot.timer < state.pot.N_rapid
                    % Case B: Rapid Sequential Saccade. 
                    % The previous saccade ended, but the eye moved again before we 
                    % could take a full stable average.
                
                    % Strategy: Look backward from now to get a common voltage point
                %    that serves as the End of the old move and Start of the new one.
                    idx_end_rapid = curr_global_idx - 1;
                    idx_start_rapid = idx_end_rapid - round(state.pot.N_avg / 2);
                    
                    mask = (state.pot.hist_index >= idx_start_rapid) & (state.pot.hist_index <= idx_end_rapid);
                    vals = state.pot.hist_values(mask);
                    
                    if ~isempty(vals)
                        V_rapid_common = mean(vals);
                        
                        % Calculate change for the previous interrupted saccade
                        if isnan(state.pot.last_peak_level)
                            diff_val = V_rapid_common - state.pot.pending_V_initial;
                        else
                            diff_val = V_rapid_common - state.pot.last_peak_level;
                        end
                       
                        current_V_est = current_V_est + diff_val;
                        
                        % Update Reference: This common point is now the "stable" level
                        state.pot.last_peak_level = V_rapid_common;
                        
                        % For the next saccade starting right now, we temporarily
                        % track V_initial, but our "Smart Logic" above will likely 
                        % ignore it in favor of last_peak_level anyway.
                        state.pot.pending_V_initial = V_rapid_common;
                    end
                    state.pot.is_waiting = false;
                else
                    % Case C: Timer was high enough, treat as normal start.
                    state.pot.is_waiting = false;
                end
            end
            
            if ~state.pot.is_waiting
                % Calculate V_initial (Baseline voltage before movement starts).
                % Note: If Smart Logic applies (2nd saccade onwards), this value 
                % is actually ignored in favor of 'last_peak_level', but we calculate 
                % it for the very first saccade or resets.
                idx_end   = curr_global_idx - 1;
                idx_start = idx_end - state.pot.N_avg;
                
                mask = (state.pot.hist_index >= idx_start) & (state.pot.hist_index <= idx_end);
                vals = state.pot.hist_values(mask);
                if ~isempty(vals)
                    state.pot.pending_V_initial = mean(vals);
                else
                    state.pot.pending_V_initial = output.V_LP(i);
                end
            end
        end
        
        % -----------------------------------------------------------
        % 3. SACCADE END (-1)
        % -----------------------------------------------------------
        if edge_val == -1
            % Movement stopped. Enter "Waiting" mode to let the signal 
            % settle (glissade).
            state.pot.is_waiting = true;
            state.pot.pending_end_idx = curr_global_idx;
            state.pot.timer = 0;
        end
    end
    % Save the updated estimate back to state and output
    state.pot.V_est = current_V_est;
    output.V_est = state.pot.V_est;
    
    %% Task 6: Convert the potentials to gaze angles
    theta = state.calib.grad * state.pot.V_est;   % [deg], from cumulative V_est
    
    
    %% Task 7: Auto recalibration (dead-zone around 0°)
    % If the estimated angle is within +/- theta_min degrees, we treat it as 0°
    % and also adjust V_est so that the internal state is re-centered (drift correction).
    
    if state.recalib.active && abs(theta) < state.recalib.theta_min
        % Remove the small drift by shifting V_est so the corresponding angle becomes 0
        % theta = grad * V_est  ->  V_est = theta / grad
        state.pot.V_est = state.pot.V_est - theta / state.calib.grad;
        theta = 0;   % snapped-to-zero angle
    end
    
    % --- Update output structure (nach der Korrektur!) ---
    output.V_est   = state.pot.V_est;
    output.ang_est = theta;
    

elseif strcmpi(cmd,'init')
    %% Task 3: Initialize the filter coefficients & filter state variables for EOG filtering
    fs = state.fs;

    % Design main LOWPASS filter (Chebyshev Type II, 35 Hz)
    fc_lp = 35; % Cutoff frequency in Hz
    Wn_lp = fc_lp / (fs/2); % Normalize to Nyquist
    [b_LP, a_LP] = cheby2(4, 60, Wn_lp, 'low');

    % Initialize filter state
    zi_LP = zeros(max(length(b_LP), length(a_LP))-1, 1);

    % Design BASELINE removal filter (Butterworth, 0.16 Hz)
    fc_baseline = 0.16; % Very low cutoff
    Wn_baseline = fc_baseline / (fs/2);
    [b_baseline, a_baseline] = butter(2, Wn_baseline, 'low');

    % Initialize filter state
    zi_baseline = zeros(max(length(b_baseline), length(a_baseline))-1, 1);
    
    %% Task 4: Parameters, Coefficients and Save Variables for Saccade Edge Detection
    % Answer:
        % With conservative big threshold, the dots are moving to the center of the step increase 
        % Too small resulted in the opposite and many wrong microsaccades were detected
        % The longer the moving average is the more points were selected to
        % calculate the average pf the derivatives
        % Skip buffers is used to skip oscillation from the beggining of the recording. If the value was too high, 
        % the saccades were detected late and the whole processing was wrong. The opposite when the value was too high.
    % Optimized for Clara Test
    state.saccade.ma_length   = 10;         % moving average length
    state.saccade.deriv_thresh = 0.005;     % derivative threshold
    state.saccade.skip_buffers = 50;         % initial buffers to ignore

    state.saccade.ma_coeffs = ones(1, state.saccade.ma_length) / state.saccade.ma_length;
    state.saccade.filter_state = zeros(state.saccade.ma_length-1, 1); % column
    state.saccade.prev_sample  = 0;
    state.saccade.prev_binary  = 0;

    %% Task 5: Parameters
    % Answer: 
    % We see a baseline dirft of an error propagated through the
    % system by the differential computation of the estimted potentials
    % The artifact free estimate is much closer to the actual eye movements
    % and does not involve microsaccades or rapid movements while
    % maintaining important information
    % The T values introduce a latency to the system. Since we must wait 
    % for T sloppy and potentially T avg to pass before confirming a 
    % saccade, the beamformer will always lag behind the actual eye 
    % movement by at least that duration (~200–300 ms).
    % Based on the course manual, putting out a single value of the 
    % artefact-free estimate per function call has the following effect on 
    % the effective sample rate: It drastically reduces the effective sample 
    % rate, effectively downsampling the signal to the rate at which the 
    % buffers are processed.

    % Timing constants
    T_avg    = 0.100; 
    T_sloppy = 0.200; 
    T_rapid  = 0.30;
    
    state.pot.N_avg    = round(T_avg * fs);
    state.pot.N_sloppy = round(T_sloppy * fs);
    state.pot.N_rapid  = round(T_rapid * fs);
    
    % Buffer setup
    state.pot.max_hist = state.pot.N_sloppy + state.pot.N_avg + fs*4; 
    state.pot.hist_values = [];
    state.pot.hist_index  = [];
    
    % State variables
    state.pot.V_est = 0;             
    state.pot.is_waiting = false;     
    state.pot.timer = 0;              
    state.pot.pending_end_idx = 0;   
    
    % NEW: Track the Last Stable Peak to avoid settling errors
    state.pot.last_peak_level = NaN; % NaN indicates we haven't seen a saccade yet
    state.pot.pending_V_initial = 0;
    
    %% Task 6: Calibration gradient
    % We already loaded calib_grad in the framework and stored it in state.calib_grad
    if isfield(state, 'calib_grad')
        state.calib.grad = state.calib_grad;
        fprintf('Using calibration gradient from framework: %g deg/unit\n', state.calib.grad);
    else
        state.calib.grad = 1;  % neutral fallback
        warning('No calib_grad found in state. Using grad = 1 (angles = potentials).');
    end

    fprintf('Using calibration gradient: %g deg/unit\n', state.calib.grad);

    %% Task 7: Parameters for auto-recalibration
    % If the estimated gaze angle is within +/- theta_min degrees,
    % we assume the user is looking at 0° and we re-center the estimate.
    
    % If the gaze angle estimate is close to 0°, we assume the user is looking straight ahead
    % and we re-center V_est. This cancels accumulated drift over time. The recalibration
    % only changes the internal offset, not the saccade detection itself.
    state.recalib.active    = true;   % Toggle on/off
    state.recalib.theta_min = 5;      % [deg] dead-zone around 0°
    
    %% Initialization Biopac (Do NOT alter)
    % The calibration gradient (calib_grad) converts the artefact-free EOG potential into gaze angle
    % The gradient is obtained from the calibration sequence and is subject-specific.
    % When plotting the angle together with the EOG signal, we may multiply it by a scale factor (×2, ×10, etc.), but
    % this does NOT change the algorithm or the actual angle estimation, only the visualization.
    % The scaling is purely to make the plotted angle curve comparable in range to the EOG signal.

    % Multiplying the angle by 10 makes its variations clearly distinguishable from the EOG potential estimate in the same plot. 
    % Without this scaling, the angle curve appears too small and visually overlaps with the EOG signal, making it difficult to see the 
    % step changes.
   
    if state.is_online
        dirpath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2 Education'; % path of mpdev header file
        biopacAPI(state.is_online,'initMPDevCom',dirpath); % initialize libraries
    else       
        biopacAPI(state.is_online,'initMPDevCom',state.EOG_file);% initialize libraries
    end
    
    biopacAPI(state.is_online, 'help');                     % print available dll functions
    biopacAPI(state.is_online, 'connectMPDev');             % Connect with MP unit
    biopacAPI(state.is_online, 'setSampleRate', state.fs);  % Set sampling rate to 500 Hz
    biopacAPI(state.is_online, 'setAcqChannels',[1 0 0 0]); % Set acquisition channels
    biopacAPI(state.is_online, 'startMPAcqDaemon');         % Start acquisition daemon
    % --> MP device is now ready to record
    biopacAPI(state.is_online, 'startAcquisition');         % called in the end to reduce delaytime upon first buffer pull
    
    output = [];
else
    error('cmd not recognized')
end
