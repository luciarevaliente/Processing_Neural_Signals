classdef dB_panning < audioPlugin
    % DB_PANNING Plugin for Chapter 3 Task 9 [cite: 581]
    % Implements a real-time panning algorithm with gain control in dB.
    
    properties
        % The parameter controlled via the slider in the GUI [cite: 581]
        % Initial value is set to 0 dB
        Gain = 0; 
    end
    
    properties (Access = private)
        % Private properties can be stored here (not required for this specific task)
        % [cite: 607]
    end
    
    methods
        % The main processing function called for each buffer [cite: 615]
        function out = process(plugin, in)
            % 'in' is the input audio buffer matrix: [Samples x Channels]
            % This task assumes stereo input (2 channels).
            
            % Initialize output buffer with the same size as input
            out = zeros(size(in));
            
            % Implementation of the panning logic from Task 9, Step 4[cite: 597]:
            % "for negative slider values, the left channel is increased in amplitude... 
            % while the right channel is decreased"
            % "vice-versa for positive slider values"
            
            % Determine gain in dB for each channel based on the slider 'Gain'
            % If Gain is -5: Left becomes +5 dB, Right becomes -5 dB.
            % If Gain is +5: Left becomes -5 dB, Right becomes +5 dB.
            gainLeft_dB  = -plugin.Gain; 
            gainRight_dB =  plugin.Gain;
            
            % Convert dB to linear scaling factor using formula: factor = 10^(dB/20)
            % [cite: 585]
            linFactorLeft  = 10^(gainLeft_dB / 20);
            linFactorRight = 10^(gainRight_dB / 20);
            
            % Apply the calculated gains to the respective channels
            out(:, 1) = in(:, 1) * linFactorLeft;  % Left Channel (Column 1)
            out(:, 2) = in(:, 2) * linFactorRight; % Right Channel (Column 2)
        end
    end
    
    % Interface methods required by the audioPlugin class [cite: 630-663]
    methods
        function out = getparamnames(plugin)
            % Returns the label for the slider in the GUI [cite: 650]
            out = {'Panning Gain (dB)'};
        end
        
        function out = getparamranges(plugin)
            % Sets the slider range from -5 to +5 as requested in Task 9, Step 1 [cite: 582]
            out = {[-5 5]};
        end
        
        function out = getnuminchan(plugin)
            % Define number of input channels (Stereo) [cite: 630]
            out = 2;
        end
        
        function out = getnumoutchan(plugin)
            % Define number of output channels (Stereo) [cite: 636]
            out = 2;
        end
    end
end