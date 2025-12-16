classdef DSB < audioPlugin
    % TASK 12 - Two-microphone Delay-and-Sum Beamformer:
        % SteeringAngle in degrees via slider (-90..+90)
        % Causal: delays only (no "look-ahead")
        % Integer-sample delay (rounded); simple and robust for the lab
        % Outputs mono beamformer result to BOTH headphone channels

    properties
        % TO SET!!!
        SteeringAngle = 0;   % degrees, slider-controlled
        MicSpacing = 0.08;   % meters 
    end

    properties (Access = private)
        fs = 44100;          % sampling rate (Hz)
        c  = 344;            % speed of sound (m/s), approx.
        maxDelaySamples = 0; % maximum delay in samples
        delayBuffer          % [maxDelaySamples x 2], stores previous samples
    end

    methods
        function initialize(plugin, initdata)
            % Determine sample rate if provided by framework; else default
            if nargin >= 2 && isstruct(initdata) && isfield(initdata, 'SampleRate') ...
                    && ~isempty(initdata.SampleRate)
                plugin.fs = initdata.SampleRate;
            else
                plugin.fs = 44100;
            end

            % Max possible inter-mic delay magnitude at |sin(phi)|=1:
            % tau_max = d/c  [seconds]
            tau_max = abs(plugin.MicSpacing) / plugin.c;

            % Convert to samples (add a couple of samples as safety margin)
            plugin.maxDelaySamples = max(1, ceil(tau_max * plugin.fs) + 2);

            % Initialize delay buffer (history)
            plugin.delayBuffer = zeros(plugin.maxDelaySamples, 2);
        end

        function out = process(plugin, in, param)
            % in: [N x 2] microphone buffer
            % out: [N x 2] stereo playback buffer (same signal in both ears)

            if size(in,2) ~= 2
                error('DSB expects exactly 2 input channels (two microphones).');
            end
            N = size(in,1);

            % Get steering angle (framework may pass param or use property)
            if nargin >= 3 && ~isempty(param)
                phi_deg = param(1);
            else
                phi_deg = plugin.SteeringAngle;
            end
            phi = deg2rad(max(-90, min(90, phi_deg)));

            % For 2 mics, relative delay magnitude:
            % tau = (d/c) * sin(phi)
            % We'll implement causally by delaying ONE channel by |tau|
            tau = (plugin.MicSpacing / plugin.c) * sin(phi);
            D = round(abs(tau) * plugin.fs); % integer delay in samples
            D = min(D, plugin.maxDelaySamples); % clamp for safety

            % Build extended buffer so we can index into "past" samples
            ext = [plugin.delayBuffer; in]; % [maxDelaySamples + N x 2]

            % Apply delay causally:
                % If tau >= 0: delay channel 2 by D, channel 1 no delay
                % If tau <  0: delay channel 1 by D, channel 2 no delay
            if tau >= 0
                x1 = ext(plugin.maxDelaySamples + 1 : plugin.maxDelaySamples + N, 1);
                x2 = ext(plugin.maxDelaySamples + 1 - D : plugin.maxDelaySamples + N - D, 2);
            else
                x1 = ext(plugin.maxDelaySamples + 1 - D : plugin.maxDelaySamples + N - D, 1);
                x2 = ext(plugin.maxDelaySamples + 1 : plugin.maxDelaySamples + N, 2);
            end

            % Delay-and-sum with gains A = 1/Nch = 1/2 (Task 12 suggestion)
            y = 0.5 * x1 + 0.5 * x2;

            % Optional safety: avoid clipping if input is hot
            y = max(min(y, 1), -1);

            % Play to both headphone channels
            out = [y, y];

            % Update delay buffer for next call
            plugin.delayBuffer = ext(end - plugin.maxDelaySamples + 1 : end, :);
        end
    end

    % Required plugin interface methods
    methods
        function out = getparamnames(plugin)
            out = {'Steering angle (deg)'};
        end

        function out = getparamranges(plugin)
            out = {[-90 90]};
        end

        function out = getnuminchan(plugin)
            out = 2;
        end

        function out = getnumoutchan(plugin, ~)
            out = 2;
        end
    end
end
