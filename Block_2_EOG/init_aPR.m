function [aPR] = init_aPR(BufferSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script initialize the sound card to be used with the Matlab Audio
% Processing Toolbox (Matlab R2017a).
%
% Obtains the list of audio devices and initializes the C400 device
% with 2 playback and 1 recording channels at a sampling rate of 44.1 kHz.
%
% Input parameters:
%       BufferSize: 512[default]
% 
% Output parameter:
%       audioPlayerRecorder object: aPR
%
% This function was written for the course Signal Processing for Audio
% Technology. Author: N. Kolotzek, Audio Information Processing, TUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create an audioPlayerRecorder object to simultaneously play and record

% default buffer size: 512 samples
if nargin<1, BufferSize=512; end

% name = 'Fast Track C400 ASIO (64-bit)';
if ispc
    %name_aPR = 'Fast Track C400 ASIO (64-bit)';
    name_aPR='Focusrite USB ASIO';
elseif ismac
    name_aPR = 'JackAudio';
else
    name_aPR = 'Fast Track C400: USB Audio (hw:1,0)';
end

aPR = audioPlayerRecorder('Device',name_aPR,...
    'SampleRate',44100,...
    'BitDepth','16-bit integer',...
    'SupportVariableSize',true,...
    'BufferSize',BufferSize,...
    'PlayerChannelMapping',[1,2],...
    'RecorderChannelMapping',1);

disp(['Setting up ' name_aPR char(10) ...
    'SampleRate: ' num2str(aPR.SampleRate) char(10) ...
    'BitDepth: ' aPR.BitDepth char(10) ...
    'BufferSize: ' num2str(aPR.BufferSize) char(10) ...
    'PlayerChannelMapping: [', num2str(aPR.PlayerChannelMapping) ']' char(10) ...
    'RecorderChannelMapping: ', num2str(aPR.RecorderChannelMapping)]);
end