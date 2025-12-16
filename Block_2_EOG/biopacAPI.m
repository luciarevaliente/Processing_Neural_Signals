function [retval ,varargout] = biopacAPI(online_offline_Flag, command, varargin)
% BIOPACAPI(online_offline_Flag, command, varargin) is an interface for 
% communication between Biopac MP unit and MATLAB. It also includes an MP
% unit simulation for offline processing. In case of unexpected behaviour 
% or error occurance biopacAPI disconnects the MP unit.
% Standard procedure:
% To obtain samples first initialize the communication, then connect the MP
% device. When sample rate and acquisition channels are set you can
% initialize an acquisition daemon which performs sampling in the
% background. Start the acquisition before pulling the first samples. After
% aquisition is finished stop the acquisition routine and cleanly
% disconnect the MP connection.
%
% Output:
%   retval     -> returns "MPSUCCESS" in case of success and an error
%                 label otherwise
%   varargout  -> as specified by "command"
%
% Input:
%   online_offline_Flag: 
%       1      -> communicate with hardware
%       2      -> offline simulation
%   command: 
%       "help"             -> show help text for specified topic
%                             varargin={} -> shows topics
%                             varargin={commandStr} -> detailed help for specified "command"
%       "initMPDevCom"     -> initializes MP communication by loading
%                             the communication library
%                             varargin={path to mpdev.h header file} -> Online Processing
%                             varargin={file name and path to EOG recording} -> Offline Processing
%       "connectMPDev"     -> connects MATLAB to the MP unit. varargin={
%                             mptype (default:103),mpmethod (default:10),
%                             serialNumber (default:"auto")}
%                             varargout={serialNumber}
%       "setSampleRate"    -> sets the sampling rate for all channels.
%                             varargin={sampleRate}.
%       "setAcqChannels"   -> define sampling channels. varargin={4
%                             elements vector, where 1 activates the
%                             channel with corresponding id}
%       "startMPAcqDaemon" -> set up a continous sampling routine running
%                             on the MP unit (with defined channels and
%                             sampling rate)
%       "startAcquisition" -> start continuous sampling on MP unit
%       "receiveMPData"    -> fetch sample buffer of all channels. call
%                             only if data acquisition has been started.
%                             varargin={number_of_samples_to_read_in}
%       "stopAcquisition"  -> stops the acquisition daemon
%       "disconnectMPDev"  -> disconnects the MP unit and terminates the
%                             communication
%
% A deeper insight in how to use this API can be gained using the "help" command.
%
% Author: Nicolas Scheiner - 15.08.2016

    %% help command for detailed function overview
    if strcmp(command, 'help')
        if nargin < 3  % show command overview
            fprintf('For more detailed help on any function enter corresponding name, e.g. biopacAPI(0,''help'',''initMPDevCom'').\nAvailable Commands:\n');
            fprintf('''initMPDevCom''\n');
            fprintf('''connectMPDev''\n');
            fprintf('''setSampleRate''\n');
            fprintf('''setAcqChannels''\n');
            fprintf('''startMPAcqDaemon''\n');
            fprintf('''startAcquisition''\n');
            fprintf('''receiveMPData''\n');
            fprintf('''stopAcquisition''\n');
            fprintf('''disconnectMPDev''\n');
        elseif strcmp(varargin{1}, 'initMPDevCom')
            fprintf('''initMPDevCom'' initializes MP communitcation by loading the communication library (online mode)\nor initializing a simulated MP device to communicate with (offline mode).\n');
            fprintf('input:\n');
            fprintf('    varargin{1}: path to ''mpdev.h'' --> eg. ''C:\Program Files (x86)\Biopac'' (ONLINE MODE)\n');
            fprintf('    varargin{1}: file location for EOG recording --> eg. ''path\eog_recording.mat'' (OFFLINE MODE)\n');
            fprintf('                 the file should contain the variable ''data'' with the EOG samples as column vectors for each channel\n');
        elseif strcmp(varargin{1}, 'connectMPDev')
            fprintf('''connectMPDev'' connects MP unit to computer.\n');
            fprintf('input:\n');
            fprintf('    varargin{1}: mptype (MP150=101, MP35=102, MP36=103->default)\n');
            fprintf('    varargin{2}: mpmethod (USB=10->default, UDP=11)\n');
            fprintf('    varargin{3}: serial number of MP150 ("auto" finds first one -> default)\n');
            fprintf('output:\n');
            fprintf('    varargout{1} : serial number of device\n');
        elseif strcmp(varargin{1}, 'setSampleRate')
            fprintf('''setSampleRate'' sets the sampling rate for all MP channels.\n');
            fprintf('input:\n');
            fprintf('    varargin{1}: sampleRate [1/sec]\n\n');
            fprintf('IMPORTANT: MP36 supports the following sampling rates ONLY\n\n');
            fprintf('Sample Rates\n');
            fprintf('[Hz]   |  [msec/sample]\n');
            fprintf('-------|----------------\n');
            fprintf('100K   |  1000/100000.0\n');
            fprintf('50K    |  1000/50000.0\n');
            fprintf('25K    |  1000/25000.0\n');
            fprintf('20K    |  1000/20000.0\n');
            fprintf('10K    |  1000/10000.0\n');
            fprintf('5000.0 |  1000/5000.0\n');
            fprintf('2000.0 |  1000/2000.0\n');
            fprintf('1000.0 |  1000/1000.0\n');
            fprintf('500.0  |  1000/500.0\n');
            fprintf('200.0  |  1000/200.0\n');
            fprintf('100.0  |  1000/100.0\n');
            fprintf('50.0   |  1000/50.0\n');
            fprintf('20.0   |  1000/20.0\n');
            fprintf('10.0   |  1000/10.0\n');
            fprintf('5.0    |  1000/5.0\n');
            fprintf('2.0    |  1000/2.0\n');
            fprintf('1.0    |  1000/1.0\n');       
        elseif strcmp(varargin{1}, 'setAcqChannels')
            fprintf('''setAcqChannels'' set acquisition channels of MP unit\n'); 
            fprintf('input:\n');
            fprintf('    varargin{1}: size(1,4) vector containing ones for active channels and zero otherwise\n'); 
        elseif strcmp(varargin{1}, 'startMPAcqDaemon')
            fprintf('''startMPAcqDaemon'' sets up a continous sampling routine running on the MP unit.\n');
            fprintf('Requires defined channels and sampling rate in order to function.\n');
        elseif strcmp(varargin{1}, 'startAcquisition')
            fprintf('''startAcquisition'' starts continuous sampling via acquisition daemon.\n');
        elseif strcmp(varargin{1}, 'receiveMPData')
            fprintf('''receiveMPData'' fetches a sample buffer of all acquisition channels. Call requires started acquisition daemon.\n');
            fprintf('input:\n');
            fprintf('    varargin{1}: amount_of_samples_to_read_in\n');  
            fprintf('output:\n');
            fprintf('    varargout{1} : buffer containing defined number of samples [mV].\n'); 
            fprintf('                   the samples are presented in an interleaved way, e.g. for channel 1 and 3 active:\n');
            fprintf('                   output = [s1_ch1, s1_ch3, s2_ch1, s2_ch3, s3_ch1 ...]\n');
        elseif strcmp(varargin{1}, 'stopAcquisition')
            fprintf('''stopAcquisition'' stops the acquisition daemon.\n');
        elseif strcmp(varargin{1}, 'disconnectMPDev')
            fprintf('''disconnectMPDev'' disconnects the MP unit, terminates the communication to MP unit, all sampling routines and deletes all previously made settings.\n');
        end
        retval = 'MPSUCCESS';   % always return 'success' for help mode

    %% online eog processing
    elseif online_offline_Flag  % online eog processing

        libname = 'mpdev'; % libname to call
        
        % initializes MP communitcation by loading library
        if strcmp(command, 'initMPDevCom')
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            end
            dll = strcat(varargin{1},'\x64\mpdev.dll');     % dll path
            doth = strcat(varargin{1},'\mpdev.h');          % h filename
            if exist(dll) ~= 3 && exist(dll) ~= 2; error('DLL file does not exist'); end;   % check for file existance
            if exist(doth) ~= 2; error('Header file does not exist'); end;                  % check for file existance
            %check if the library is already loaded
            if libisloaded(libname)
                fprintf('Library still loaded from old session. Restoring default state...\n');
                biopacAPI(online_offline_Flag, libname, 'disconnectMPDev');
                unloadlibrary(libname);
            end
            warning off MATLAB:loadlibrary:enumexists;      % turn off annoying enum warnings
            loadlibrary(dll,doth);                          % load the library
            fprintf('MP communication initialized. mpdev.dll loaded\n');                    % indicate successful function behaviour
            
        % connect MP unit to computer
        % input:
        %   varargin{1}: mptype (MP150=101, MP35=102, MP36=103->default)
        %   varargin{2}: mpmethod (USB=10->default, UDP=11)
        %   varargin{3}: serial number of MP150 ("auto" finds first one -> default)
        % output:
        %   varargout{1} : serial number of device
        elseif strcmp(command, 'connectMPDev')
            % check input
            if nargin < 3
                varargin{1} = 103;
                varargin{2} = 10;
                varargin{3} = 'auto';
            elseif nargin < 4
                varargin{2} = 10;
                varargin{3} = 'auto';
            elseif nargin < 5
                varargin{3} = 'auto';
            end
            [retval, varargout{1}] = calllib(libname,'connectMPDev',varargin{1},varargin{2},varargin{3});
            if testMPError(online_offline_Flag, retval, 'Failed to Connect.'); return; end;
            fprintf('MP unit connected\n');           % indicate successful function behaviour
            
        % set the sampling rate for all MP channels
        % input:
        %   varargin{1}: sampleRate [1/sec]
        %
        % IMPORTANT: MP36 supports the following sampling rates ONLY
        %
        % Sample Rates
        % [Hz]   |  [msec/sample]
        % -------|----------------
        % 100K	 |  1000/100000.0
        % 50K	 |  1000/50000.0
        % 25K	 |  1000/25000.0
        % 20K	 |  1000/20000.0
        % 10K	 |  1000/10000.0
        % 5000.0 |	1000/5000.0
        % 2000.0 |	1000/2000.0
        % 1000.0 |	1000/1000.0
        % 500.0	 |  1000/500.0
        % 200.0	 |  1000/200.0
        % 100.0	 |  1000/100.0
        % 50.0	 |  1000/50.0
        % 20.0	 |  1000/20.0
        % 10.0	 |  1000/10.0
        % 5.0	 |  1000/5.0
        % 2.0	 |  1000/2.0
        % 1.0	 |  1000/1.0
        elseif strcmp(command, 'setSampleRate')
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            % valid sampling freq?
            elseif (varargin{1} ~= 1     && varargin{1} ~= 2     && varargin{1} ~= 5    && ...
                    varargin{1} ~= 10    && varargin{1} ~= 20    && varargin{1} ~= 50   && ...
                    varargin{1} ~= 100   && varargin{1} ~= 200   && varargin{1} ~= 500  && ...
                    varargin{1} ~= 1000  && varargin{1} ~= 2000  && varargin{1} ~= 5000 && ...
                    varargin{1} ~= 10000 && varargin{1} ~= 20000 && varargin{1} ~= 50000 && varargin{1} ~= 100000)
                if testMPError(online_offline_Flag, 'MPERROR', 'Invalid sampling rate. Proper rates can be found in the documentation.'); retval = 'MPERROR'; return; end;
            end
            retval = calllib(libname, 'setSampleRate', 1000/varargin{1});  % varargin{1} = sample duration in msec
            if testMPError(online_offline_Flag, retval, 'Failed to Set Sample Rate.'); return; end;
            fprintf('Sample Rate set to %d Hz\n',varargin{1});      % indicate successful function behaviour
        
        % set acquisition channels of MP unit
        % input:
        %   varargin{1}: size(1,4) vector containing 1 for active channels and zero otherwise
        elseif strcmp(command, 'setAcqChannels')
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            end            
            retval = calllib(libname, 'setAcqChannels',int32(varargin{1})); % dll takes int32 4x1 vector (varargin(id) = 1 indicates channel id is active)
            if testMPError(online_offline_Flag, retval, 'Failed to Set Acq Channels.'); return; end;
            fprintf('Aquisition Channels set to: '); fprintf('Ch.%d ',find(varargin{1})); fprintf('\n'); % indicate successful function behaviour
        
        % set up a continous sampling routine running on the MP unit.
        % works only with defined channels and sampling rate
        elseif strcmp(command, 'startMPAcqDaemon')
            retval = calllib(libname, 'startMPAcqDaemon');
            if testMPError(online_offline_Flag, retval, 'Failed to Start Acquisition Daemon.'); return; end;
            fprintf('Acquisition Daemon started\n');              % indicate successful function behaviour
    
        % starts continuous sampling via acquisition daemon
        elseif strcmp(command, 'startAcquisition')
            retval = calllib(libname, 'startAcquisition');
            if testMPError(online_offline_Flag, retval, 'Failed to Start Acquisition.'); return; end;
            fprintf('Acquisition started\n');                     % indicate successful function behaviour
        
        % fetch a sample buffer of all acquisition channels. Call requires started acquisition daemon.
        % input:
        %     varargin{1}: amount_of_samples_to_read_in
        % output:
        %     varargout{1} : buffer containing defined number of samples.
        %                    the samples are presented in an interleaved way, e.g. for channel 1 and 3 active:
        %                    output = [s1_ch1, s1_ch3, s2_ch1, s2_ch3, s3_ch1 ...]
        elseif strcmp(command, 'receiveMPData')
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            end  
            [retval, varargout{1}]  = calllib(libname, 'receiveMPData',zeros(1,varargin{1}), varargin{1}, 0); 
            varargout{1} = varargout{1} * 1000; % conversion to mV
            if testMPError(online_offline_Flag, retval, 'Failed to receive MP data.'); return; end;
        
        % stop the acquisition daemon
        elseif strcmp(command, 'stopAcquisition')
            retval = calllib(libname, 'stopAcquisition');
            if testMPError(online_offline_Flag, retval, 'Failed to Stop.'); return; end;
            fprintf('Acquisition stopped\n');                     % indicate successful function behaviour
        
        % disconnect the MP unit
        elseif strcmp(command, 'disconnectMPDev')
            retval = calllib(libname, 'disconnectMPDev');
            if testMPError(online_offline_Flag, retval, 'Failed to Disconnect.'); return; end;
            if libisloaded(libname); unloadlibrary(libname); end; % unload library
            fprintf('MP unit disconnected\n');                    % indicate successful function behaviour
                
        else % if command cannot be found
            retval = 'MPERROR';
            if testMPError(online_offline_Flag, retval, 'MP Command not found.'); return; end;
            
        end
        
    %% offline eog processing
    else    % offline eog processing
        
        persistent MPDev; % simulated MPDev object - persistent across calls to this function
            
        % initializes MP communitcation by loading eog recording and
        % setting up simulated MP device
        if strcmp(command, 'initMPDevCom')
            if ~exist(varargin{1},'file'); error('EOG mat-file does not exist'); end; % check for file existance
            load(varargin{1}); % --> eg. 'path\some_eog_recording.mat' the file should contain the variable
                               %     "data" with the EOG samples as column vectors for each channel
            % check file content
            if ~exist('data','var') % check existance
                error('EOG mat-file does not contain a variable called "data"');
            elseif size(data,2) > size(data,1) % check size
                data = data'; % transpose data in case it in the wrong format
            end
            MPDev = struct('connectionFlag',0,'sampleRate',0, ...
                'acquChannels',[0 0 0 0],'daemonFlag',0, ...
                'daemonStartedFlag',0,'bufferOffset',0, ...
                'time',[],'data',data); % create simulated MP unit
            fprintf('MP communication initialized\n');                  % indicate successful function behaviour
            
        % connect MP unit to computer
        % input:
        %   varargin{1}: mptype (MP150=101, MP35=102, MP36=103->default)
        %   varargin{2}: mpmethod (USB=10->default, UDP=11)
        %   varargin{3}: serial number of MP150 ("auto" finds first one -> default)
        % output:
        %   varargout{1} : serial number of device
        elseif strcmp(command, 'connectMPDev')
            try 
                MPDev.connectionFlag = MPDev.connectionFlag + 1; % set connection = 1; -> create error if access not possible
            catch
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires initialization first.']); retval = 'MPERROR'; return; end;
            end
            fprintf('MP unit connected\n');          % indicate successful function behaviour
            
        % set the sampling rate for all MP channels
        % input:
        %   varargin{1}: sampleRate [1/sec]
        %
        % IMPORTANT: MP36 supports the following sampling rates ONLY
        %
        % Sample Rates
        % [Hz]   |  [msec/sample]
        % -------|----------------
        % 100K	 |  1000/100000.0
        % 50K	 |  1000/50000.0
        % 25K	 |  1000/25000.0
        % 20K	 |  1000/20000.0
        % 10K	 |  1000/10000.0
        % 5000.0 |	1000/5000.0
        % 2000.0 |	1000/2000.0
        % 1000.0 |	1000/1000.0
        % 500.0	 |  1000/500.0
        % 200.0	 |  1000/200.0
        % 100.0	 |  1000/100.0
        % 50.0	 |  1000/50.0
        % 20.0	 |  1000/20.0
        % 10.0	 |  1000/10.0
        % 5.0	 |  1000/5.0
        % 2.0	 |  1000/2.0
        % 1.0	 |  1000/1.0            
        elseif strcmp(command, 'setSampleRate')
            % check working connection
            if ~MPDev.connectionFlag
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires connection first.']); retval = 'MPERROR'; return; end;
            end
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            % valid sampling freq?
            elseif (varargin{1} ~= 1     && varargin{1} ~= 2     && varargin{1} ~= 5    && ...
                    varargin{1} ~= 10    && varargin{1} ~= 20    && varargin{1} ~= 50   && ...
                    varargin{1} ~= 100   && varargin{1} ~= 200   && varargin{1} ~= 500  && ...
                    varargin{1} ~= 1000  && varargin{1} ~= 2000  && varargin{1} ~= 5000 && ...
                    varargin{1} ~= 10000 && varargin{1} ~= 20000 && varargin{1} ~= 50000 && varargin{1} ~= 100000)
                if testMPError(online_offline_Flag, 'MPERROR', 'Invalid sampling rate. Proper rates can be found in the documentation.'); retval = 'MPERROR'; return; end;
            end
            MPDev.sampleRate = varargin{1}; % set sample rate of simulated MP unit
            fprintf('Sample Rate set to %d Hz\n',varargin{1});      % indicate successful function behaviour
        
        % set acquisition channels of MP unit
        % input:
        %   varargin{1}: size(1,4) vector containing 1 for active channels and zero otherwise    
        elseif strcmp(command, 'setAcqChannels')
            % check working connection
            if ~MPDev.connectionFlag
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires connection first.']); retval = 'MPERROR'; return; end;
            end
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            end            
            MPDev.acquChannels = varargin{1}; % set acquisition channels of simulated MP unit
            fprintf('Aquisition Channels set to: '); fprintf('Ch.%d ',find(varargin{1})); fprintf('\n'); % indicate successful function behaviour
        
        % set up a continous sampling routine running on the MP unit -> set state
        % works only with defined channels and sampling rate
        elseif strcmp(command, 'startMPAcqDaemon')
            % check state of simulated MP unit
            if ~MPDev.connectionFlag || ~find(MPDev.acquChannels) || ~MPDev.sampleRate
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires connection and set parameters for acquisition channels and parameters first.']); retval = 'MPERROR'; return; end;
            end
            MPDev.daemonFlag = 1;                                   % set AcqDaemon flag of simulated MP unit
            fprintf('Acquisition Daemon started\n');                % indicate successful function behaviour
        
        % starts continuous sampling via acquisition daemon -> timer setup
        elseif strcmp(command, 'startAcquisition')
            % check state of simulated MP unit
            if ~MPDev.connectionFlag || ~MPDev.daemonFlag
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires an acquisition daemon first.']); retval = 'MPERROR'; return; end;
            end
            MPDev.daemonStartedFlag = 1;                            % set AcqDaemonStarted flag of simulated MP unit
            MPDev.time = tic;                                       % start timer for simulated acquisition routine
            fprintf('Acquisition started\n');                       % indicate successful function behaviour
        
        % fetch a sample buffer of all acquisition channels. Call requires started acquisition daemon.
        % input:
        %     varargin{1}: amount_of_samples_to_read_in
        % output:
        %     varargout{1} : buffer containing defined number of samples.
        %                    the samples are presented in an interleaved way, e.g. for channel 1 and 3 active:
        %                    output = [s1_ch1, s1_ch3, s2_ch1, s2_ch3, s3_ch1 ...]
        elseif strcmp(command, 'receiveMPData')
            % check state of simulated MP unit
            if ~MPDev.daemonStartedFlag
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires a started acquisition daemon.']); retval = 'MPERROR'; return; end;
            end
            
            time = toc(MPDev.time);                                     % time since acquisition started
            availableSample = time*MPDev.sampleRate-MPDev.bufferOffset; % calc amount of samples in buffer
            if availableSample < varargin{1}    % if amount of samples in buffer < samples requested
                pause((varargin{1}-availableSample)/MPDev.sampleRate);  % wait until available
                fprintf('Waiting for buffer...\n');                     % indicate timing problem
            end
            % check input
            if nargin < 3
                if testMPError(online_offline_Flag, 'MPERROR', ['Not enough arguements. ', command, ' requires 1 argument.']); retval = 'MPERROR'; return; end;
            end            
                    %old version with 2 input arguments: load(libname); % --> eg. 'some_eog_recording.mat'
                    %                                    varargout{1} = data((varargin{2}-1)*varargin{1}+1:varargin{1}*varargin{2},1)'; % return buffer #varargout{2} with size varargout{1}
            try
                varargout{1} = MPDev.data(MPDev.bufferOffset+1:MPDev.bufferOffset+varargin{1},1)'; % return buffer with size varargout{1}
            catch
%                 if testMPError(online_offline_Flag, 'MPERROR', 'Recording not long enough. Execution stopped.'); retval = 'MPERROR'; return; end;
                  % End of recording reached, ouptut zeros
                    varargout{1} = zeros(1,varargin{1});  
            end
            MPDev.bufferOffset = MPDev.bufferOffset + varargin{1}; % update buffer offset of simulated MP unit
        
        % stop the acquisition daemon
        elseif strcmp(command, 'stopAcquisition')
            if MPDev.daemonStartedFlag
                MPDev.daemonStartedFlag = 0;                        % reset AcqDaemonStarted flag of simulated MP unit
            else
                if testMPError(online_offline_Flag, 'MPERROR', [command, ' requires a started acquisition daemon.']); retval = 'MPERROR'; return; end;
            end
            fprintf('Acquisition stopped\n');                       % indicate successful function behaviour
            
        % disconnect the MP unit, terminates all sampling routines and deletes all previously made settings    
        elseif strcmp(command, 'disconnectMPDev')
            if MPDev.connectionFlag
                MPDev.conncetionFlag = 0;                           % reset simulated MP unit
                MPDev.sampleRate = 0;                               % reset simulated MP unit
                MPDev.acquChannels = [0 0 0 0];                     % reset simulated MP unit
                MPDev.bufferOffset = 0;                           	% reset simulated MP unit
                MPDev.daemonFlag = 0;                               % reset simulated MP unit
                MPDev.daemonStartedFlag = 0;                        % reset simulated MP unit
                MPDev.time = [];                                    % reset simulated MP unit
                MPDev.data = [];                                    % reset simulated MP unit
            end
            fprintf('MP unit disconnected\n');                      % indicate successful function behaviour
            
        else % if command cannot be found
            retval = 'MPERROR';
            if testMPError(online_offline_Flag, retval, 'MP Command not found.'); retval = 'MPERROR'; return; end;
            
        end
        
        retval = 'MPSUCCESS';   % always return 'success' for offline mode
        
    end
    
end

%% test function -> confirm correct libcall
% TESTMPERROR is called after each call to the biopac dll to confirm the
% correct command execution. Returns int(0) for correct execution and
% int(1) and an error message otherwise
function errFlag = testMPError(online_offline_Flag, retval, errMsg)
   if ~strcmp(retval,'MPSUCCESS')
       errFlag = 1; % return error
       fprintf('%s\n',errMsg); % display error message
       if online_offline_Flag 
           calllib('mpdev', 'disconnectMPDev');
       end
   else
       errFlag = 0; % return "no error"
   end
end