% clear the workspace
close all;
clear all
clc;

%% General Intialization Parameters (ADJUST VALUES HERE ONLY) %%%%%%%%%%%%%%
% fixed variables
fs_acq = 500;               % sampling frequency [Hz]
Nacq_buff = 6;              % EOG acquisition buffer size

% variables to be set
runtime = 25;               % set here the runtime of the processing loop in seconds
plot_flag = 1;              % plot on = 1; off = 0;
plot_freq = 1;              % Frequency of plotting [# of acquisition buffers] --> NECESSARY TO BE 1 IN ORDER TO PLOT ALL THE SCATTER DOTS
is_online = 0;              % online = 1; offline = 0;
pause_time = 0.02;  %0.02

% Paths for Task 3
% state.EOG_file = 'C:\Users\User\Enginyeria de dades - UAB\4t TUM\Neural_signal_lab\Block_2_Repo\Lucia_Table_1_2.mat';
% state.EOG_file = 'C:\Users\User\Enginyeria de dades - UAB\4t TUM\Neural_signal_lab\Block_2_Repo\Lucia_Table_2.mat';
state.EOG_file = 'EOGTestClara.mat';

% Fixing the offset of the offline EOG file (only when offline)
if ~is_online
    raw = load(state.EOG_file);
    signal = raw.data;

    % Remove DC offset
    offset        = (max(signal) + min(signal)) / 2;
    signal_zeroed = signal - offset;

    % Normalize to +1/-1
    signal_fixed  = signal_zeroed / max(abs(signal_zeroed));

    % Save corrected signal into a new file as variable "data"
    data = signal_fixed;  % keep the same variable name "data" for compatibility
    [p,n,e] = fileparts(state.EOG_file);
    corrected_file = fullfile(p, [n '_corrected' e]);
    save(corrected_file, 'data');

    % Tell the rest of the framework to use the corrected file
    state.EOG_file = corrected_file;
end


% Path for Task 6
% calib_file     = 'C:\Users\User\Enginyeria de dades - UAB\4t TUM\Neural_signal_lab\Block_2_Repo\Lucia_Calibration.mat';
calib_file = 'C:\Users\User\Enginyeria de dades - UAB\4t TUM\Neural_signal_lab\Block_2_Repo\grad_calib.mat';

calib_data        = load(calib_file); % Load calibration gradient (variable "calib_grad")
state.calib_grad  = calib_data.calib_grad;   % g=10


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize RT_EOG

% write the variables needed for initialization to state
state.fs = fs_acq;                
state.Nacq_buff = Nacq_buff;              
state.is_online = is_online;

[~,state] = RT_EOG('init',state);

% total number of samples
nTotal = ceil(state.fs*runtime/state.Nacq_buff)*state.Nacq_buff;

% buffers for all the signal
EOG_Vraw = zeros(nTotal,1);
EOG_VLP  = zeros(nTotal,1);
EOG_VBP  = zeros(nTotal,1);
EOG_Vest = zeros(nTotal,1);   % artefact-free cumulative EOG estimate
EOG_ang  = zeros(nTotal,1);   % eye-gaze angle estimate [deg]

% Plotting Routine 1/5: setting up the figure===============================
if plot_flag
    figure
    xlabel('Time (s)')
    set(gca,'Xlim',[0, runtime],'Ylim',[-3.5 5]) %runtime
    hold on
end
%===========================================================================
%% Real-Time Processing Loop
% here you define how many buffers you fetch 
for fr_idx = 1:ceil(state.fs*runtime/state.Nacq_buff)
    state.fr_idx = fr_idx;
    [output,state] = RT_EOG('process',state);    

    idx = (fr_idx-1)*state.Nacq_buff + (1:state.Nacq_buff);

    % save all
    EOG_Vraw(idx) = output.V_raw;
    EOG_VLP(idx)  = output.V_LP;
    EOG_VBP(idx)  = output.V_BP;
    EOG_Vest(idx) = output.V_est;  % scalar estimate per buffer → broadcast to all samples in this buffer
    EOG_ang(idx)  = ones(state.Nacq_buff, 1) * output.ang_est;
    
    %% Plotting Routine 2/5 (Task 3) ======================================================
    % if plot_flag && ~mod(fr_idx,plot_freq)
    %     idx_offset = (fr_idx-plot_freq)*state.Nacq_buff;
    %     plot_range = idx_offset + (1:plot_freq*state.Nacq_buff);
    % 
    %     t_range = plot_range/state.fs;
    % 
    %     % choose what to see in real-time
    %     % plot(t_range, EOG_Vraw(plot_range),'b')   % cruda
    %     % plot(t_range, EOG_VLP(plot_range),'r')     % low-pass
    %     plot(t_range, EOG_VBP(plot_range),'k')    % band-pass
    % end


    %% Plotting Routine 3/5 (Task 4)======================================================
    % if plot_flag && ~mod(fr_idx,plot_freq)
    %     idx_offset = (fr_idx-plot_freq)*state.Nacq_buff;
    %     plot_range = idx_offset + (1:plot_freq*state.Nacq_buff);
    %     t_range    = plot_range / state.fs;
    % 
    %     % --- Plot Low-Pass EOG ---
    %     plot(t_range, EOG_VLP(plot_range), 'r'); hold on;
    % 
    %     % --- Extract edge indices from THIS buffer ---
    %     edges_range = output.edge_idx;                  % vector of length Nacq_buff
    %     edges_idx   = find(edges_range ~= 0);           % positions inside buffer
    %     edges_type  = edges_range(edges_idx);           % +1 or -1
    % 
    %     % --- Convert buffer-local indices to global sample indices ---
    %     global_edges = (fr_idx-1)*state.Nacq_buff + edges_idx;
    %     t_edges      = global_edges / state.fs;
    % 
    %     % --- Plot start (+1) and end (-1) markers ---
    %     if any(edges_type == 1)
    %         scatter(t_edges(edges_type==1), ...
    %              EOG_VLP(global_edges(edges_type==1)),  "green");
    %              % 'g*', 'MarkerSize', 1, 'LineWidth', 2); % START (green)
    %     end
    % 
    %     if any(edges_type == -1)
    %         %fprintf("End of Saccade\n");
    %         scatter(t_edges(edges_type==-1), ...
    %              EOG_VLP(global_edges(edges_type==-1)),  "blue");
    %              % 'r*', 'MarkerSize', 1, 'LineWidth', 2); % END (red)
    %     end
    % 
    % end


    %% Plotting Routine 4/5 (Task 5) - Add V_est ======================================================
    % if plot_flag && ~mod(fr_idx,plot_freq)
    %     idx_offset = (fr_idx-plot_freq)*state.Nacq_buff;
    %     plot_range = idx_offset + (1:plot_freq*state.Nacq_buff);
    %     t_range = plot_range / state.fs;
    % 
    %     % Plot Low-Pass EOG
    %     plot(t_range, EOG_VLP(plot_range), 'r'); hold on;
    % 
    %     % --- CHANGED SECTION START ---
    %     % Calculate current time
    %     t_current = fr_idx * state.Nacq_buff / state.fs;
    % 
    %     % Initialize plotting state variables if they don't exist yet
    %     % (We store these in 'state' so they persist between loop iterations)
    %     if ~isfield(state, 'plot_last_t')
    %          state.plot_last_t = max(0, t_current - (state.Nacq_buff*plot_freq/state.fs));
    %          state.plot_last_V = 0; 
    %     end
    % 
    %     % Plot connected Step Function:
    %     % 1. Horizontal line from Last Time to Current Time (holding previous value)
    %     % 2. Vertical line up/down to the New Value
    %     plot([state.plot_last_t, t_current, t_current], ...
    %          [state.plot_last_V, state.plot_last_V, output.V_est], ...
    %          'b-', 'LineWidth', 1); % Thinner line
    % 
    %     % Update plotting history for the next loop
    %     state.plot_last_t = t_current;
    %     state.plot_last_V = output.V_est;
    %     % --- CHANGED SECTION END ---
    % 
    %     % Plot saccade markers (green/blue dots)
    %     edges_range = output.edge_idx;
    %     edges_idx = find(edges_range ~= 0);
    %     edges_type = edges_range(edges_idx);
    %     global_edges = (fr_idx-1)*state.Nacq_buff + edges_idx;
    %     t_edges = global_edges / state.fs;
    %     if any(edges_type == 1)
    %         scatter(t_edges(edges_type==1), ...
    %              EOG_VLP(global_edges(edges_type==1)), "green");
    %     end
    %     if any(edges_type == -1)
    %         scatter(t_edges(edges_type==-1), ...
    %              EOG_VLP(global_edges(edges_type==-1)), "blue");
    %     end
    % end
    
    
    %% Plotting Routine 5/5 (Task 6) - LP + V_est + angle ======================================================
    if plot_flag && ~mod(fr_idx, plot_freq)
    
        idx_offset = (fr_idx-plot_freq)*state.Nacq_buff;
        plot_range = idx_offset + (1:plot_freq*state.Nacq_buff);
        t_range    = plot_range / state.fs;
    
        % 1) Low-pass EOG
        plot(t_range, EOG_VLP(plot_range), 'r'); 
        hold on;
    
        % Current buffer time
        t_current = fr_idx * state.Nacq_buff / state.fs;
    
        % 2) Initialize plotting state once
        if ~isfield(state, 'plot_last_t')
            state.plot_last_t   = 0;
            state.plot_last_V   = output.V_est;
            state.plot_last_ang = output.ang_est;
        end
    
        % --- SCALE ANGLE VISUALLY BY 10 ---
        % ×10 scaling improves interpretability, not accuracy. The algorithm always uses the true (unscaled) angle; the extra gain 
        % is only for clearer plotting.
        angle_plot_prev = state.plot_last_ang * 10;
        angle_plot_curr = output.ang_est * 10;
        % -----------------------------------
    
        % 3) Step for V_est (blue)
        plot([state.plot_last_t, t_current, t_current], ...
             [state.plot_last_V, state.plot_last_V, output.V_est], ...
             'b-', 'LineWidth', 1);
    
        % 4) Step for ANGLE (magenta, 10× scaled)
        plot([state.plot_last_t, t_current, t_current], ...
             [angle_plot_prev, angle_plot_prev, angle_plot_curr], ...
             'm--', 'LineWidth', 1);
    
        % 5) Update memory
        state.plot_last_t   = t_current;
        state.plot_last_V   = output.V_est;
        state.plot_last_ang = output.ang_est;
    
        % 6) Saccade markers
        edges_range = output.edge_idx;
        edges_idx   = find(edges_range ~= 0);
        edges_type  = edges_range(edges_idx);
    
        global_edges = (fr_idx-1)*state.Nacq_buff + edges_idx;
        t_edges      = global_edges / state.fs;
    
        if any(edges_type == 1)
            scatter(t_edges(edges_type==1), ...
                EOG_VLP(global_edges(edges_type==1)), "green");
        end
        if any(edges_type == -1)
            scatter(t_edges(edges_type==-1), ...
                EOG_VLP(global_edges(edges_type==-1)), "blue");
        end
    
    end

    pause(pause_time);
end


%% Stop BIOPAC
biopacAPI(state.is_online,'stopAcquisition')
biopacAPI(state.is_online,'disconnectMPDev')


%% Final Step: Save Results
% Store filtered signals and gaze angle estimates after all tasks are finished

save_filename = ['EOG_Results_' datestr(now, 'yyyymmdd_HHMMSS') '.mat'];

fprintf('Saving processed data to %s...\n', save_filename);

save(save_filename, ...
    'EOG_Vraw', ...   % Raw EOG signal
    'EOG_VLP',  ...   % Low-pass filtered EOG
    'EOG_VBP',  ...   % Band-pass filtered EOG
    'EOG_Vest', ...   % Artefact-free potential estimate
    'EOG_ang', ...    % Final gaze angle estimate
    'fs_acq', ...
    'Nacq_buff');