function varargout = FeedBackLoopBoLi_20190622(varargin)

% 20220708 update: add keyboard interrupt to feedback loop
% 20220606 update: in feedback loop main, change plotting from AWG
% electronic in to AWG optical out
% 20220527 update: add movie capture function in feedback loop main
% 20220524 update: AWG on range changed to 1-6 V. Not working yet. 
% 
% modified by SZhao on 12062020 to initially try to adapt to 114 MHz. 
% 
% FEEDBACKLOOPBOLI_20190622 MATLAB code for FeedBackLoopBoLi_20190622.fig
%      FEEDBACKLOOPBOLI_20190622, by itself, creates a new FEEDBACKLOOPBOLI_20190622 or raises the existing
%      singleton*.
%
%      H = FEEDBACKLOOPBOLI_20190622 returns the handle to a new FEEDBACKLOOPBOLI_20190622 or the handle to
%      the existing singleton*.
%
%      FEEDBACKLOOPBOLI_20190622('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEEDBACKLOOPBOLI_20190622.M with the given input arguments.
%]

%      FEEDBACKLOOPBOLI_20190622('Property','Value',...) creates -a new FEEDBACKLOOPBOLI_20190622 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FeedBackLoopBoLi_20190622_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FeedBackLoopBoLi_20190622_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FeedBackLoopBoLi_20190622

% Last Modified by GUIDE v2.5 21-Jun-2019 15:07:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FeedBackLoopBoLi_20190622_OpeningFcn, ...
                   'gui_OutputFcn',  @FeedBackLoopBoLi_20190622_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% 0.0 Parameters define
% --- Executes just before FeedBackLoopBoLi_20190622 is made visible.
function FeedBackLoopBoLi_20190622_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
% varargin   command line arguments to FeedBackLoopBoLi_20190622 (see VARARGIN)

% Choose default command line output for FeedBackLoopBoLi_20190622
handles.output = hObject;

% Initialization status
set(handles.edit21, 'String', 'Waiting for initialization'); 

% AWG
% default AWG ID
set(handles.edit20, 'String', 'Slot5'); 
% default gain, 0 to 6 V
set(handles.edit6, 'String', num2str(1)); 
% default offset, -0.25 to 0.25 V
set(handles.edit7, 'String', num2str(-0.015)); 
% default sampling rate: 0 - 114 MHz
set(handles.edit15, 'String', num2str(114)); 
% default sampling source: 1. Internal 2. External
set(handles.text26, 'String', 'External'); 
% default sampling source: 1. Internal 2. External
set(handles.edit43, 'String', num2str(0)); 
% default repetition rate, MHz
set(handles.edit46, 'String', num2str(21)); 


% OSC
% default average times
set(handles.edit28, 'String', num2str(1));   % in MHz
% default sampling clock source
set(handles.text134, 'String', 'External');   % 1.Internal 2. External
% default sampling rate, 1.Internal: 500 or 1000 MHz, 2. External: 912 MHz
set(handles.edit11, 'String', num2str(912));   % in MHz
% default number of buffers in one acquisition 
set(handles.edit12, 'String', num2str(1));      % no unit
% default buffer duration
set(handles.edit16, 'String', num2str(701));   % in us
% default board ID
set(handles.edit22, 'String', num2str(1));      % no unit
% default trigger level
set(handles.edit27, 'String', num2str(100));      % in mV, -400 to 400 mV
% default trigger source, it can be 1. Channel A, 2. Channel B, 3. External .
set(handles.text44, 'String', 'External '); 
% Show result or not? 1. Yes, 2. No .(No+space)
set(handles.text48, 'String', 'No '); 
% Available channels? 1. A  , 2. B  ,3. A+B.
set(handles.text50, 'String', 'A  '); 
% boardHandle
addpath('AlazarDriver'); 
AlazarDefs; 
if ~alazarLoadLibrary() 
    fprintf('Error: ATSApi library not loaded\n');
    return
end
systemId = int32(1); % 1 is ID for the device
boardId = int32(1);  % 1 is ID for the device
boardHandle = AlazarGetBoardBySystemID(systemId, boardId);
setdatatype(boardHandle, 'voidPtr', 1, 1);
handles.boardHandle = boardHandle;


% Feedback
% default loop number
set(handles.edit17, 'String', num2str(1)); 
% default work mode: 1. Test mode 2. Work mode
% set(handles.text24, 'String', 'Work mode');
% Channel A starts at ... us    When see spikes on the right edge of PD
% envelopes, this number needs to decrease
% set(handles.edit29, 'String', '44.75'); % for 4 MHz
set(handles.edit29, 'String', '0.17'); % for 32 MHz
% Channel B starts at ... us
% set(handles.edit30, 'String', '42.86'); % for 4 MHz
set(handles.edit30, 'String', '42.88'); % for 32 MHz
% Pause for each loop
set(handles.edit31, 'String', '0.001');
% Image to pattern, Yes or No .
set(handles.text94, 'String', 'No ');
% show raw data, Yes or No .
set(handles.text108, 'String', 'No ');
% Pulse peak test
set(handles.edit53, 'String', '0');
% Waveform show x limit, left side
set(handles.edit58, 'String', '0');
% Waveform show x limit, right side
set(handles.edit56, 'String', '139.9');
% Waveform show x offset
set(handles.edit57, 'String', '1');
% Figure? What figure is going to show? 1. Simple, 2. Full  . 3. No    .
set(handles.text144, 'String', 'Full  ');


% Data Save
% default Date: it is the current day
set(handles.edit1, 'String', num2str(yyyymmdd(datetime)));
% default Group
set(handles.edit2, 'String', 'Test group');
% default Object
set(handles.edit3, 'String', '1');
% default folder, pwd for current folder
set(handles.edit23, 'String', pwd);
% default if save data:
handles.DataSaveSaveModeNum = 1;
NewNum = rem(handles.DataSaveSaveModeNum,2)+1;
handles.DataSaveSaveMode = ['Yes';'No '];
set(handles.text25, 'String', handles.DataSaveSaveMode(NewNum,:));
set(handles.text36, 'String', 'Not save yet');
set(handles.text37, 'String', 'Not save yet');
set(handles.text38, 'String', 'Not save yet');
% default Folder for read
set(handles.edit33, 'String', 'F:\Bo\Resonant_scanner\Aribitrary_pattern_20180322');
% default Folder for read
set(handles.edit34, 'String', 'OSC_Measure_20180123_Data1_1');


% Image process
handles.ImageProcessFolder = 'F:\NI_PXI5422_GUI_20190621\Pattern_generation';
set(handles.edit25, 'String', handles.ImageProcessFolder);
handles.ImageProcessFileName = '66kHz_522by522_strip1.tif';
set(handles.edit24, 'String', handles.ImageProcessFileName);
handles.captionFontSize = 8;
% default Threshold
set(handles.edit35, 'String', '140');
% default target grows by ... pixels
set(handles.edit36, 'String', '0');
% default image filter ... pixels
set(handles.edit37, 'String', '100');
% default file type, Stack or Image
set(handles.text70, 'String', 'Stack');
% default normalization factor
set(handles.edit38, 'String', '0.7');
% default pulse value in scan back time
set(handles.edit39, 'String', '0.6');
% default compensation pixels for blank line
set(handles.edit40, 'String', '10');
% default value: the flyback time is divided by this factor
set(handles.edit41, 'String', '1');
% default pixel time
set(handles.edit42, 'String', '40'); % 32 MHz laser, 1 pixel 1 pulse
% default if equal each line, 1. Yes. 2. No .
set(handles.text104, 'String', 'Yes');
set(handles.edit47, 'String', '140'); % Scanback time, left side
set(handles.edit48, 'String', '1566'); % Scan time
set(handles.edit49, 'String', '140'); % Scan back time, right side
set(handles.edit50, 'String', '0.0'); % Scan back time, right side2
set(handles.edit51, 'String', '0'); % Scan back time extra compensation for feedback, a line has 2606 AWG 4047 samples, now it has floor((583+1441*2+560)/4)*4+4 = 4028, so compensation = 4047-4028 = 19.
set(handles.edit52, 'String', '0'); % Fly back time, % max time = 33100 samples!! reasonable values: [29000, 33100]
set(handles.edit54, 'String', '20230503_lineRate66k_90FF.txt'); % Image pixels to AWG samples.

% For 114 MHz below:
% set(handles.edit47, 'String', '583'); % Scanback time, left side
% set(handles.edit48, 'String', '1441'); % Scan time
% set(handles.edit49, 'String', '560'); % Scan back time, right side
% set(handles.edit50, 'String', '0.0'); % Scan back time, right side2
% set(handles.edit51, 'String', '19'); % Scan back time extra compensation for feedback, a line has 2606 AWG 4047 samples, now it has floor((583+1441*2+560)/4)*4+4 = 4028, so compensation = 4047-4028 = 19.
% set(handles.edit52, 'String', '32000'); % Fly back time, % max time = 33100 samples!! reasonable values: [29000, 33100]
% set(handles.edit54, 'String', 'Pixel_vs_samples20190118.txt'); % Image pixels to AWG samples.


% Initial AWG
[hObject,handles] = AWGInitializationFunction(hObject, eventdata, handles);


% AWG Initialization status update
AWGInitializationCheck(handles);
% OSC Initialization status update
OSCInitializationCheck(handles);

set(handles.edit21, 'String', 'GUI is ready  ^_^'); 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FeedBackLoopBoLi_20190622 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FeedBackLoopBoLi_20190622_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

if get(handles.text22, 'String')=='Initialized'
	message = sprintf('The AWG is not removed yet.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Remove and exit', 'Exit without removing', 'Yes');
	if strcmpi(reply, 'Remove and exit')
       [hObject, handles]=AWGCheckAndRemove(hObject, eventdata, handles); % AWG remove
    end
end
% 
% In case the GUI does not close, run the two lines below
%  set(0,'ShowHiddenHandles','on')
% delete(get(0,'Children'))
% 
% Update handles structure
guidata(hObject, handles);


delete(hObject);


% 1.1 AWG initialization check
function AWGInitializationCheck(handles)
% Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    set(handles.text22, 'String', 'Initialized'); 
else
    set(handles.text22, 'String', 'Not initialized'); 
end

% 1.2 AWG initialization
function [NewhObject,Newhandles] = AWGInitializationFunction(hObject, ~, handles)
% AWG Check and remove
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
    disconnect(handles.AWGdeviceObj);
    delete(handles.AWGdeviceObj);
    clear handles.AWGdeviceObj;
    handles = rmfield(handles,'AWGdeviceObj');
    disp('AWG is turned off and deleted.');
end

% Create a MATLAB Instrument Object
resourceID=get(handles.edit20, 'String');
handles.AWGdeviceObj = icdevice('niFgen',resourceID); % Create an instrument object from the niFgen driver
connect(handles.AWGdeviceObj);                        % Connect driver instance
% fopen(handles.AWGdeviceObj)
disp('AWG is initialized.');
%Query maximum capabilities
% [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities');
% disp('Max waveform size is ' + string(MAXIMUMWAVEFORMSIZE) + newline + 'Min size is ' + string(MINIMUMWAVEFORMSIZE));
Newhandles=handles;
NewhObject=hObject;

% 1.3 AWG check and remove
function [NewhObject, Newhandles]=AWGCheckAndRemove(hObject, ~, handles)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
    disconnect(handles.AWGdeviceObj);
    delete(handles.AWGdeviceObj);
    clear handles.AWGdeviceObj;
    handles = rmfield(handles,'AWGdeviceObj');
    disp('AWG is turned off and deleted.');
end
NewhObject=hObject;
Newhandles=handles;

% 1.4 AWG run -- modified by SZhao 12082020
function AWG8on8off(~, ~, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Status check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% AWG Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice ~= 1
    errordlg('AWG is not intialized','GUI Error');
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters define %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters define, these values are defined in the driver's header file 'niFgen.h'
NIFGEN_VAL_OUTPUT_ARB = 1;      % Arbitrary mode
NIFGEN_VAL_OUTPUT_SEQ = 2;
ChannelName = '0';              % This value is described in the help file 'NI Signal Generators Help'
Enabled = 1;                    % Enable output waveform
% Rising_Edge = 101;              % Start trigger: rising edge detection
% Falling_Edge = 102;             % Start trigger: falling edge detection
% NIFGEN_VAL_CONTINUOUS = 2;      % Trigger mode: continuous 
% NIFGEN_VAL_STEPPED = 3;         % Trigger mode: stepped 
NIFGEN_VAL_MARKER_EVENT = 1001; % Marker event, for Marker 0

% AWG group
AWGSamplingRate = str2double(get(handles.edit15, 'String'))*1e6;     % Sampling rate, E15, in Hz
    
NIFGEN_ATTR_ARB_GAIN = str2double(get(handles.edit6, 'String'));     % Gain, E6, Range (0,6)
    if NIFGEN_ATTR_ARB_GAIN>6 || NIFGEN_ATTR_ARB_GAIN<=0
        errordlg('Gain range is [0, 6]','GUI Error');
        return;
    end
NIFGEN_ATTR_ARB_OFFSET = str2double(get(handles.edit7, 'String'));   % Offset, E7, Range: (-0.25,0.25)
    if NIFGEN_ATTR_ARB_OFFSET>0.25 || NIFGEN_ATTR_ARB_OFFSET<-0.25
        errordlg('AWG offset range is [-0.25, 0.25]','GUI Error');
        return;
    end
SamplingSource = get(handles.text26, 'String');                   % Sampling source, T26, P14, P15, 1. Internal 2. External
    if SamplingSource == 'Internal'
        if AWGSamplingRate>114e6 || AWGSamplingRate<=0.01e6
            errordlg('Sampling range is [0.01, 114] MHz','GUI Error');
            return;
        end
        SamplingSourceNum = 0;
    else
        if SamplingSource == 'External'
            SamplingSourceNum = 1;
        else
            errordlg('Sampling source of AWG is not defined','GUI Error');
            return;
        end
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% AWG waveform define %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseOn = 1;                % Waveform value when pulse is on
PulseOff = 0;               % Waveform value when pulse is off


DTime = 1/AWGSamplingRate;      % s
TWindow = DTime*32;           % s
NSample = floor(TWindow/DTime/4)*4;                    % 8 = 32 MHz / 4 MHz; One pulse, 8 samples.
% TWindow/DTime

% Pattern 1: 10 neurons (200 um) in one line(500 um)
% TWindowOn = TWindow*50/100;                            % Time window for pulse on
% NSampleOn = floor(TWindowOn/DTime/8)*8;                 % Sample number for pulse on
% WAVEFORMOn = linspace(PulseOn,PulseOn,NSampleOn);       % Waveform for pulse on
% NSampleOff = NSample-NSampleOn;                         % Sample number for pulse off
% WAVEFORMOff = linspace(PulseOff,PulseOff,NSampleOff);   % Waveform for pulse off
% WAVEFORMDATAARRAY = [WAVEFORMOn,WAVEFORMOff];          % Waveform, value: (-1,1)
% % WAVEFORMDATAARRAYNorm = abs((WAVEFORMDATAARRAY-PulseOff)/max(WAVEFORMDATAARRAY-PulseOff)); % Normalized
% WAVEFORMSIZE = length(WAVEFORMDATAARRAY);               % Waveform size
% % WaveformIdealOut = WAVEFORMDATAARRAY;                   % Ideal final optical output waveform, normalized.
% % WaveformIdealMatrix = reshape(WAVEFORMDATAARRAY,[1,WAVEFORMSIZE/1]); % Matrix
% % WaveformIdealPeak = abs((WaveformIdealMatrix-PulseOff)/max(WaveformIdealMatrix-PulseOff));           % Ideal final output peak.
% % NumPulseOnOne = sum(WaveformIdealPeak);                    % number of pulse on
% WAVEFORMDATAARRAY2 = [WAVEFORMOn,WAVEFORMOff];          % Waveform, value: (-1,1)
% WAVEFORMSIZE2 = length(WAVEFORMDATAARRAY2);               % Waveform size

% if SamplingSourceNum==1 % for external sampling source
    NSample = 50;
    NSampleOn = 25;                                        % Sample number for pulse on 
    WAVEFORMOn = linspace(PulseOn,PulseOn,NSampleOn);       % Waveform for pulse on
    NSampleOff = NSample - NSampleOn;                      % Sample number for pulse off
    WAVEFORMOff = linspace(PulseOff,PulseOff,NSampleOff);   % Waveform for pulse off
    OnePeriodRef = [WAVEFORMOff,WAVEFORMOn];          % Waveform, value: (-1,1)
%     NSample = 4000;
%     NSampleOn = 55;                                        % Sample number for pulse on 
    WAVEFORMOn = linspace(PulseOn,PulseOn,NSampleOn);       % Waveform for pulse on
    NSampleOff = NSample - NSampleOn;                    % Sample number for pulse off
    WAVEFORMOff = linspace(PulseOff,PulseOff,NSampleOff);   % Waveform for pulse off
    OnePeriod = [WAVEFORMOff,WAVEFORMOn]; 
%     OnePeriod(NSampleOff-5409:NSampleOff-5400) = PulseOn;
    WaveformRef = repmat(OnePeriodRef,1,20);
    Waveform2nd = repmat(OnePeriod,1,20);
    WAVEFORMDATAARRAY = [WaveformRef,Waveform2nd];
    
%     WAVEFORMDATAARRAY = [WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff]; 
    % WAVEFORMDATAARRAYNorm = abs((WAVEFORMDATAARRAY-PulseOff)/max(WAVEFORMDATAARRAY-PulseOff)); % Normalized
    WAVEFORMSIZE = length(WAVEFORMDATAARRAY);               % Waveform size
    
    
%     dividor = 32;
%     dividorArray = linspace(PulseOff,PulseOff, WAVEFORMSIZE);
%     dividorArray(1:dividor:end) = 1;
%     WAVEFORMDATAARRAY = WAVEFORMDATAARRAY.*dividorArray;
    
    
    
    % WaveformIdealOut = WAVEFORMDATAARRAY;                   % Ideal final optical output waveform, normalized.
    % WaveformIdealMatrix = reshape(WAVEFORMDATAARRAY,[1,WAVEFORMSIZE/1]); % Matrix
    % WaveformIdealPeak = abs((WaveformIdealMatrix-PulseOff)/max(WaveformIdealMatrix-PulseOff));           % Ideal final output peak.
    % NumPulseOnOne = sum(WaveformIdealPeak);                    % number of pulse on
%     WAVEFORMDATAARRAY2 = [WAVEFORMOn,WAVEFORMOff];          % Waveform, value: (-1,1)
%     WAVEFORMSIZE2 = length(WAVEFORMDATAARRAY2);               % Waveform size
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Enable AWG output %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disable the AWG if it is running. (To protect AWG)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
    disp('AWG is stoped since it was on.');
end


% Configure the sampling rate
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configuresamplerate',AWGSamplingRate);
if SamplingSourceNum==1 % if use external source
    invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'configuresampleclocksource',"ClkIn");
end
% Configure output mode - SEQ
% invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_SEQ);
% % Configure output mode - ARB
invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_ARB);
% Create arbitrary sequence
% [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);
% [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE2] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE2,WAVEFORMDATAARRAY2);
% WAVEFORMHANDLESARRAY = [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE2];
% LOOPCOUNTSARRAY = [1,1];
% SequenceLength = 2;
% % [NIFGEN_ATTR_ARB_SEQUENCE_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'createarbsequence',SequenceLength,WAVEFORMHANDLESARRAY,LOOPCOUNTSARRAY);
% % Advanced arb sequence
% NIFGEN_VAL_NO_MARKER=-1;
% SAMPLECOUNTSARRAY=[WAVEFORMSIZE,WAVEFORMSIZE2];
% MARKERLOCATIONARRAY=[0,NIFGEN_VAL_NO_MARKER];
% COERCEDMARKERSARRAY=[0,0];
% [COERCEDMARKERSARRAY,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'createadvancedarbsequence',SequenceLength,WAVEFORMHANDLESARRAY,LOOPCOUNTSARRAY,SAMPLECOUNTSARRAY,MARKERLOCATIONARRAY,COERCEDMARKERSARRAY);



% % Create arbitrary waveform
[NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);



% scnim
% % Query maximum capabilities
% [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities');
% fprintf(['MAXIMUMNUMBEROFWAVEFORMS: ', num2str(MAXIMUMNUMBEROFWAVEFORMS), ' \n']);
% fprintf(['MAXIMUMWAVEFORMSIZE: ', num2str(MAXIMUMWAVEFORMSIZE), ' \n']);


% % Configure arbitrary waveform
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configurearbwaveform',ChannelName,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);
% Configure arbitrary sequency
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'configurearbsequence',ChannelName,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);






% export Marker 0 to PFI 1
set(handles.AWGdeviceObj.Arbitrarywaveformarbitrarywaveformmode(1), 'Marker_Position', 0);
invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'exportsignal',NIFGEN_VAL_MARKER_EVENT,"Marker0","PFI1");

% Tigger mode :Continuouse/Stepped
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuretriggermode',ChannelName,NIFGEN_VAL_STEPPED);
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuredigitaledgestarttrigger',"PFI0",Rising_Edge);

% Continuous mode [The default mode is continuous mode]
% NIFGEN_VAL_OPERATE_CONTINUOUS=0;
% invoke(handles.AWGdeviceObj.Configuration,'configureoperationmode',ChannelName,NIFGEN_VAL_OPERATE_CONTINUOUS)

% Initiate the Waveform Generation
invoke(handles.AWGdeviceObj.Waveformcontrol,'initiategeneration');

% Enable the Output
invoke(handles.AWGdeviceObj.Configuration,'configureoutputenabled', ChannelName, Enabled);
% fprintf(['AWG takes time: ', num2str(toc(testTime)), 's. \n']);

% disp('AWG is working');


% 1.5 AWG run 2  -- modified by SZhao on 20200430
function AWG32MHzTo4MHz(~, ~, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Status check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AWG Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice ~= 1
    errordlg('AWG is not intialized','GUI Error');
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters define %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters define, these values are defined in the driver's header file 'niFgen.h'
NIFGEN_VAL_OUTPUT_ARB = 1;      % Arbitrary mode
ChannelName = '0';              % This value is described in the help file 'NI Signal Generators Help'
Enabled = 1;                    % Enable output waveform
% Rising_Edge = 101;              % Start trigger: rising edge detection
% Falling_Edge = 102;             % Start trigger: falling edge detection
% NIFGEN_VAL_CONTINUOUS = 2;      % Trigger mode: continuous 
% NIFGEN_VAL_STEPPED = 3;         % Trigger mode: stepped 
NIFGEN_VAL_MARKER_EVENT = 1001; % Marker event, for Marker 0

% AWG group
AWGSamplingRate = str2double(get(handles.edit15, 'String'))*1e6;     % Sampling rate, E15, in Hz
    
NIFGEN_ATTR_ARB_GAIN = str2double(get(handles.edit6, 'String'));     % Gain, E6, Range (0,6)
    if NIFGEN_ATTR_ARB_GAIN>6 || NIFGEN_ATTR_ARB_GAIN<=0
        errordlg('Gain range is [0, 6]','GUI Error');
        return;
    end
NIFGEN_ATTR_ARB_OFFSET = str2double(get(handles.edit7, 'String'));   % Offset, E7, Range: (-0.25,0.25)
    if NIFGEN_ATTR_ARB_OFFSET>0.25 || NIFGEN_ATTR_ARB_OFFSET<-0.25
        errordlg('AWG offset range is [-0.25, 0.25]','GUI Error');
        return;
    end
SamplingSource = get(handles.text26, 'String');                   % Sampling source, T26, P14, P15, 1. Internal 2. External
    if SamplingSource == 'Internal'
        if AWGSamplingRate>200e6 || AWGSamplingRate<=0.01e6
            errordlg('Sampling range is [0.01, 200] MHz','GUI Error');
            return;
        end
        SamplingSourceNum = 0;
    else
        if SamplingSource == 'External'
            SamplingSourceNum = 1;
        else
            errordlg('Sampling source of AWG is not defined','GUI Error');
            return;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% AWG waveform define %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseOn = 1;                % Waveform value when pulse is on
PulseOff = 0;               % Waveform value when pulse is off

RepRateDivider = str2double(get(handles.edit46, 'String')); % Read how loops
DTime = 1/AWGSamplingRate;      % s
% test. Delete below!
% if RepRate == 66
%     WAVEFORMDATAARRAYTT = linspace(PulseOn, PulseOn,50);
%     WAVEFORMDATAARRAYTT(7:50)= PulseOff;
%     WAVEFORMDATAARRAY= [WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT];
% end
% test. Delete above.
WAVEFORMDATAARRAYTT = linspace(PulseOn, PulseOn, RepRateDivider);
WAVEFORMDATAARRAYTT(2:RepRateDivider)= PulseOff;
WAVEFORMDATAARRAY= [WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT ...
    WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT];
% if RepRateDivider == 1
%     WAVEFORMDATAARRAYTT = linspace(PulseOn, PulseOn,114);
%     WAVEFORMDATAARRAY= [WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT];
% end
% if RepRateDivider == 2
%     WAVEFORMDATAARRAYTT = linspace(PulseOn, PulseOn,114);
%     WAVEFORMDATAARRAYTT(2:2:114)= PulseOff;
%     WAVEFORMDATAARRAY= [WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT];
% end
% if RepRateDivider == 6
% WAVEFORMDATAARRAYTT = linspace(PulseOn, PulseOn,6);
% WAVEFORMDATAARRAYTT(2:6)= PulseOff;
% WAVEFORMDATAARRAY= [WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT WAVEFORMDATAARRAYTT];
% end


% TWindow = DTime*8;           % s
% NSample = round(TWindow/DTime/4)*4;                    % 8 = 32 MHz / 4 MHz; One pulse, 8 samples.
% TWindow/DTime

% % Pattern 1: 10 neurons (200 um) in one line(500 um)
% TWindowOn = TWindow/8;                            % Time window for pulse on
% NSampleOn = round(TWindowOn/DTime);                 % Sample number for pulse on
% WAVEFORMOn = linspace(PulseOn,PulseOn,NSampleOn);       % Waveform for pulse on
% NSampleOff = NSample-NSampleOn;                         % Sample number for pulse off
% WAVEFORMOff = linspace(PulseOff,PulseOff,NSampleOff);   % Waveform for pulse off
% WAVEFORMDATAARRAY = [WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff,WAVEFORMOn,WAVEFORMOff];          % Waveform, value: (-1,1)
% WAVEFORMDATAARRAYNorm = abs((WAVEFORMDATAARRAY-PulseOff)/max(WAVEFORMDATAARRAY-PulseOff)); % Normalized
WAVEFORMSIZE = length(WAVEFORMDATAARRAY);               % Waveform size
% WaveformIdealOut = WAVEFORMDATAARRAY;                   % Ideal final optical output waveform, normalized.
% WaveformIdealMatrix = reshape(WAVEFORMDATAARRAY,[1,WAVEFORMSIZE/1]); % Matrix
% WaveformIdealPeak = abs((WaveformIdealMatrix-PulseOff)/max(WaveformIdealMatrix-PulseOff));           % Ideal final output peak.
% NumPulseOnOne = sum(WaveformIdealPeak);                    % number of pulse on
% figure
% plot([WAVEFORMOn,WAVEFORMOff]);
for k=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Enable AWG output %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disable the AWG if it is running. (To protect AWG)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
    disp('AWG is stoped since it was on.');
end

% Configure the sampling rate
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configuresamplerate',AWGSamplingRate);    % Tell the AWG an expected sampling rate or the internal sampling rate
if SamplingSourceNum==1 % if use external source
    invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'configuresampleclocksource',"ClkIn");
end
% AWG output delay, this value should be less than 1 sample time
AWGADJUSTMENTTIME = str2double(get(handles.edit43, 'String'))*1e-6; % Read how loops


% Configure AWG sample delay
if AWGADJUSTMENTTIME~=0
    invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'adjustsampleclockrelativedelay',AWGADJUSTMENTTIME);
end
% set(handles.AWGdeviceObj.Clocksadvanced(1), 'External_Clock_Delay_Binary_Value', 0e-9);
% set(handles.AWGdeviceObj.Output(1), 'Channel_Delay', 0);
set(handles.AWGdeviceObj.Clocksadvanced(1), 'Sample_Clock_Absolute_Delay', 0e-9);
% Create arbitrary waveform
[NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);

%Query maximum capabilities
% [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities');
% disp('Max waveform size is ' + string(MAXIMUMWAVEFORMSIZE) + newline + 'Min size is ' + string(MINIMUMWAVEFORMSIZE));

% Configure output mode
invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_ARB);

% Configure arbitrary waveform
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configurearbwaveform',ChannelName,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);

% export Marker 0 to PFI 1
set(handles.AWGdeviceObj.Arbitrarywaveformarbitrarywaveformmode(1), 'Marker_Position', 0);
invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'exportsignal',NIFGEN_VAL_MARKER_EVENT,"Marker0","PFI1");

% Tigger mode :Continuouse/Stepped
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuretriggermode',ChannelName,NIFGEN_VAL_STEPPED);
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuredigitaledgestarttrigger',"PFI0",Rising_Edge);

% Continuous mode [The default mode is continuous mode]
% NIFGEN_VAL_OPERATE_CONTINUOUS=0;
% invoke(handles.AWGdeviceObj.Configuration,'configureoperationmode',ChannelName,NIFGEN_VAL_OPERATE_CONTINUOUS)

% Initiate the Waveform Generation
invoke(handles.AWGdeviceObj.Waveformcontrol,'initiategeneration');

% Enable the Output
invoke(handles.AWGdeviceObj.Configuration,'configureoutputenabled', ChannelName, Enabled);
% disp('AWG is working');
end

% 1.6 AWG disable
function [NewhObject, Newhandles]=AWGDisable(hObject, handles)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
end
NewhObject=hObject;
Newhandles=handles;


% 2.1 OSC initialization check
function OSCInitializationCheck(handles)

InitialStatus = 1;
% Add path to AlazarTech mfiles
% Bo: add driver path, driver is placed in the same file
addpath('AlazarDriver'); 

% Call mfile with library definitions
% Bo: definition of various parameters, this m file is in driver folder
AlazarDefs; 

% Load driver library
% Bo: if ~(m file): if m file exist, run it, if not, fprintf.
if ~alazarLoadLibrary() 
    fprintf('Error: ATSApi library not loaded\n');
    InitialStatus = 0;
end
boardHandle = handles.boardHandle;
if boardHandle.Value == 0
    fprintf('Error: Unable to open board system ID %u board ID %u\n', systemId, boardId);
    InitialStatus = 0;
end
% Configure the board's sample rate, input, and trigger settings
if ~configureBoard(handles)
    fprintf('Error: Board configuration failed\n');
    InitialStatus = 0;
end

if InitialStatus == 1
    set(handles.text23, 'String', 'Initialized'); 
else
    set(handles.text23, 'String', 'Not initialized'); 
end



% 2.2 OSC Measure
function [NewhObject, Newhandles]=OSCMeasureTest(hObject, ~, handles)
fprintf('Measure begins...\n');
% Add path to AlazarTech mfiles
% Bo: add driver path, driver is placed in the same file
% addpath('AlazarDriver'); 
% 
% % Call mfile with library definitions
% % Bo: definition of various parameters, this m file is in driver folder
% AlazarDefs; 
% 
% % Load driver library
% % Bo: if ~(m file): if m file exist, run it, if not, fprintf.
% alazarLoadLibrary();

% testTime = tic;
%%% OSC group
% Sampling rate set, edit 11, in Hz
OSCSamplingRate = str2double(get(handles.edit11, 'String'))*1e6;
if OSCSamplingRate~=500e6 && OSCSamplingRate~=1000e6 && OSCSamplingRate~= 912e6 && OSCSamplingRate~= 917e6
    errordlg('OSC sampling rate should be 500 or 1000 or 912 MHz','GUI Error');
    return;
end
% Buffer time, edit 16, in second
OSCBufferTime = str2double(get(handles.edit16, 'String'))*1e-6;
% Buffer numbers in 1 acquisition
OSCBufferNumberInOneAcq = str2double(get(handles.edit12, 'String'));
% Acquisition Time, in second
OSCAcquisitionTime = OSCBufferTime * OSCBufferNumberInOneAcq;
% Samples in one total acquisition
OSCNlengthAllChannel = OSCAcquisitionTime*OSCSamplingRate;
% Channels
OSCChannelNumber = 1;
OSCChannelNumberCode = get(handles.text50, 'String');
if OSCChannelNumberCode == 'A+B'
    OSCChannelNumber = 2;
end


% Samples in one channel in one acquisition, last 3 samples are deleted.
OSCNlengthPerChannel = floor(OSCNlengthAllChannel)-3;
% OSC time
DTimeOSC = 1/OSCSamplingRate;
VTimeOSC = linspace(1,OSCNlengthPerChannel,OSCNlengthPerChannel)*DTimeOSC;


% Configure the board's sample rate, input, and trigger settings
if ~configureBoard(handles)
    fprintf('Error: Board configuration failed\n');
    return
end

% fprintf(['OSC configuration takes time: ',num2str(toc(testTime)),' s.\n']);

% Acquire data, optionally saving it to a file
AverageTime = str2double(get(handles.edit28, 'String'));
VItensityAmp = VTimeOSC*0;
VItensityMod = VItensityAmp;



for AverageNum=1:AverageTime
    [result,OSCDataReceived]=acquireData(handles);
    if result ~=1
        fprintf('Error: Acquisition failed\n');
        return
    end
    if OSCChannelNumber == 2
        VItensityAmp = VItensityAmp+OSCDataReceived(1:2:(OSCNlengthPerChannel-1)*2+1); % Channel A, amplified signal
        VItensityMod = VItensityMod+OSCDataReceived(2:2:(OSCNlengthPerChannel)*2); % Channel B, modulated signal
    else
        VItensityAmp = VItensityAmp+OSCDataReceived(1:OSCNlengthPerChannel); % Channel A, amplified signal
    end
end

% Fourier transform

% DTimeOSC = VTimeOSC(2)-VTimeOSC(1);
% NSampOSC = length(VTimeOSC);
% FWindowOSC = 1/DTimeOSC;                         % F : N*df, frequency window
% DFreqOSC = FWindowOSC/(NSampOSC-1);                 % df : minimum frequency
% VFreqOSC = (-(NSampOSC-1)/2:(NSampOSC-1)/2)'*DFreqOSC; % f : frequency vector
% 
% VItensityAmpSpec = DTimeOSC*fftshift(fft(ifftshift(VItensityAmp)));
% FilterOSC=exp(-VFreqOSC.^2/2/(100e6)^2)';
% VItensityAmpSpecFiltered = VItensityAmpSpec.*FilterOSC;
% 
% figure
% plot(VFreqOSC/1e6,VItensityAmpSpec,'b');
% hold on
% plot(VFreqOSC/1e6,VItensityAmpSpecFiltered,'r');
% hold on
% xlim([-100,100])

% fprintf(['OSC takes time after acquisition: ',num2str(toc(testTime)),' s.\n']);

if OSCChannelNumber == 2
    VItensityAmp = VItensityAmp/AverageTime;
    VItensityMod = VItensityMod/AverageTime;
    figure
    subplot(2,1,1);
    plot(VTimeOSC/1e-6,VItensityMod);
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    title('Channel B');

    subplot(2,1,2);
    plot(VTimeOSC/1e-6,VItensityAmp);
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    % ylim([-0.05,1]);
    title('Channel A');    
else
    VItensityAmp = VItensityAmp/AverageTime;
%     VItensityAmpFiltered = 1/DTimeOSC*fftshift(ifft(ifftshift(VItensityAmpSpecFiltered)));
    figure
    plot(VTimeOSC/1e-6,VItensityAmp,'b');
    hold on
%     plot(VTimeOSC/1e-6,VItensityAmpFiltered,'r','linewidth',2);
%     hold on
%     OSCNlengthPerChannelInt = floor(OSCNlengthPerChannel/32);
%     plot(VTimeOSC(1:32:32*OSCNlengthPerChannelInt)/1e-6,VItensityAmp(1:32:32*OSCNlengthPerChannelInt),'o');
%     hold on
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    % ylim([-0.05,1]);
    title('Channel A');    
end
fprintf('Measure finished.\n');

NewhObject=hObject;
Newhandles=handles;


% 2.3 OSC Measure and save
function [NewhObject, Newhandles]=OSCMeasureAndSave(hObject, ~, handles)
fprintf('Measure begins...\n');


% Add path to AlazarTech mfiles
% Bo: add driver path, driver is placed in the same file
% addpath('AlazarDriver'); 
% 
% % Call mfile with library definitions
% % Bo: definition of various parameters, this m file is in driver folder
% AlazarDefs; 
% 
% % Load driver library
% % Bo: if ~(m file): if m file exist, run it, if not, fprintf.
% alazarLoadLibrary();

%%% OSC group
% Sampling rate set, edit 11, in Hz
OSCSamplingRate = str2double(get(handles.edit11, 'String'))*1e6;
if OSCSamplingRate~=500e6 && OSCSamplingRate~=1000e6 && OSCSamplingRate~= 912e6
    errordlg('OSC sampling rate should be 500 or 1000 or 912 MHz','GUI Error');
    return;
end
% Buffer time, edit 16, in second
OSCBufferTime = str2double(get(handles.edit16, 'String'))*1e-6;
% Buffer numbers in 1 acquisition
OSCBufferNumberInOneAcq = str2double(get(handles.edit12, 'String'));
% Acquisition Time, in second
OSCAcquisitionTime = OSCBufferTime * OSCBufferNumberInOneAcq;
% Samples in one total acquisition
OSCNlengthAllChannel = OSCAcquisitionTime*OSCSamplingRate;
% Channels
OSCChannelNumber = 1;
OSCChannelNumberCode = get(handles.text50, 'String');
if OSCChannelNumberCode == 'A+B'
    OSCChannelNumber = 2;
end
% Samples in one channel in one acquisition, last 3 samples are deleted.
OSCNlengthPerChannel = floor(OSCNlengthAllChannel)-3;
% OSC time
DTimeOSC = 1/OSCSamplingRate;
VTimeOSC = linspace(1,OSCNlengthPerChannel,OSCNlengthPerChannel)*DTimeOSC;


% Configure the board's sample rate, input, and trigger settings
if ~configureBoard(handles)
    fprintf('Error: Board configuration failed\n');
    return
end

% Acquire data, optionally saving it to a file
AverageTime = str2double(get(handles.edit28, 'String'));
VItensityAmp = VTimeOSC*0;
VItensityMod = VItensityAmp;

for AverageNum=1:AverageTime
[result,OSCDataReceived]=acquireData(handles);
if result ~=1
    fprintf('Error: Acquisition failed\n');
    return
end
if OSCChannelNumber == 2
    VItensityAmp = VItensityAmp+OSCDataReceived(1:2:(OSCNlengthPerChannel-1)*2+1); % Channel A, amplified signal
    VItensityMod = VItensityMod+OSCDataReceived(2:2:(OSCNlengthPerChannel)*2); % Channel B, modulated signal
else
    VItensityAmp = VItensityAmp+OSCDataReceived(1:OSCNlengthPerChannel); % Channel A, amplified signal
end
end

if OSCChannelNumber == 2
    VItensityAmp = VItensityAmp/AverageTime;
    VItensityMod = VItensityMod/AverageTime;
    figure
    subplot(2,1,1);
    plot(VTimeOSC/1e-6,VItensityMod);
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    title('Channel B');
 
    subplot(2,1,2);
    plot(VTimeOSC/1e-6,VItensityAmp);
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    % ylim([-0.05,1]);
    title('Channel A');  
    % Save data
    PresentFileName = ['OSC_Measure_',get(handles.edit1, 'String'),'_',get(handles.edit2, 'String'),'_',get(handles.edit3, 'String')];
    save([get(handles.edit23, 'String'),'\',PresentFileName],'VTimeOSC','VItensityMod','VItensityAmp');
    set(handles.text38, 'String',get(handles.text37, 'String'));
    set(handles.text37, 'String',get(handles.text36, 'String'));
    set(handles.text36, 'String',PresentFileName);
else
    VItensityAmp = VItensityAmp/AverageTime;
    figure
    plot(VTimeOSC/1e-6,VItensityAmp);
    xlabel('Time/us');
    ylabel('Voltage/mV');
    xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
    % ylim([-0.05,1]);
    title('Channel A');   
    % Save data
    PresentFileName = ['OSC_Measure_',get(handles.edit1, 'String'),'_',get(handles.edit2, 'String'),'_',get(handles.edit3, 'String')];
    save([get(handles.edit23, 'String'),'\',PresentFileName],'VTimeOSC','VItensityAmp');
    set(handles.text38, 'String',get(handles.text37, 'String'));
    set(handles.text37, 'String',get(handles.text36, 'String'));
    set(handles.text36, 'String',PresentFileName);
end

fprintf('Measure finished.\n');
NewhObject=hObject;
Newhandles=handles;


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, ~, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles]=FeedbackloopRun(hObject, handles);
guidata(hObject, handles);



% 3.0 Run last waveform - stepped trigger mode: using scanimage line clock
function [NewhObject, Newhandles]=RunLastWaveformStepped(hObject, handles,StepOrContinuous)
% The feedback loop function starts from here
fprintf('\n\n\n');
fprintf('Run last calculated waveform\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Status check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AWG Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice ~= 1
    errordlg('AWG is not intialized','GUI Error');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters define %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters define, these values are defined in the driver's header file 'niFgen.h'
NIFGEN_VAL_OUTPUT_ARB = 1;      % Arbitrary mode
NIFGEN_VAL_OUTPUT_SEQ = 2;
ChannelName = '0';              % This value is described in the help file 'NI Signal Generators Help'
Enabled = 1;                    % Enable output waveform
Rising_Edge = 101;              % Start trigger: rising edge detection
Falling_Edge = 102;             % Start trigger: falling edge detection
NIFGEN_VAL_CONTINUOUS = 2;      % Trigger mode: continuous 
NIFGEN_VAL_STEPPED = 3;         % Trigger mode: stepped 
NIFGEN_VAL_MARKER_EVENT = 1001; % Marker event, for Marker 0
% FirstTest =0;


%%% AWG group
% Sampling rate, E15, in Hz
AWGSamplingRate = str2double(get(handles.edit15, 'String'))*1e6;     
if AWGSamplingRate>114e6 || AWGSamplingRate<=0.01e6
    errordlg('Sampling range is [0.01, 114] MHz','GUI Error');
    return;
end
% Gain, E6, Range (0,6)
NIFGEN_ATTR_ARB_GAIN = str2double(get(handles.edit6, 'String'));
if NIFGEN_ATTR_ARB_GAIN>6 || NIFGEN_ATTR_ARB_GAIN<=0
    errordlg('Gain range is [0, 6]','GUI Error');
    return;
end
% Offset, E7, Range: (-0.25,0.25)
NIFGEN_ATTR_ARB_OFFSET = str2double(get(handles.edit7, 'String'));   
if NIFGEN_ATTR_ARB_OFFSET>0.25 || NIFGEN_ATTR_ARB_OFFSET<-0.25
    errordlg('AWG offset range is [-0.25, 0.25]','GUI Error');
    return;
end
% Sampling source, T26, P14, P15, 1. Internal 2. External
AWGSamplingSource = get(handles.text26, 'String');      
% AWG output delay, this value should be less than 1 sample time
AWGADJUSTMENTTIME = str2double(get(handles.edit43, 'String'))*1e-6; % Read how loops


WAVEFORMDATAARRAY = handles.LastWAVEFORMDATAARRAY;          % Waveform, value: (-1,1)
WAVEFORMSIZE = length(WAVEFORMDATAARRAY);

PatternInSampleFirst = handles.FinalPatternInSampleFirst ;
PatternInSampleLast  = handles.FinalPatternInSampleLast  ;
    % First part, line 2 to 511
    LineOfThePatternFirst = length(PatternInSampleFirst(:,1));
    SamplesOfOneLineFirst = length(PatternInSampleFirst(1,:));
    % Second part, line 512 to line 1
    SamplesOfOneLineLast = length(PatternInSampleLast);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Enable AWG output %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Disable the AWG if it is running. (To protect AWG)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
end



% [hObject,handles] = AWGInitializationFunction(hObject, 1, handles);
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'cleararbmemory');

% set(handles.AWGdeviceObj.Output(1), 'Output_Enabled', 0.0);
% invoke(handles.AWGdeviceObj.Waveformcontrol,'abortgeneration');
    

% if kLoopnum>1
% %     Configurationfunctionsarbitrarysequenceoutput
%     invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);
% end

% Sample rate multiplier: may be this is for output sample clock!!!
% set(handles.AWGdeviceObj.Clocksadvanced, 'External_Sample_Clock_Multiplier', 3)
% get83 = get(handles.AWGdeviceObj.Clocksadvanced, 'External_Sample_Clock_Multiplier')

% Configure the sampling rate
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configuresamplerate',AWGSamplingRate);
if AWGSamplingSource=='External' % The next line is only for external source.
   invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'configuresampleclocksource',"ClkIn");
end

% Configure AWG sample delay
if AWGADJUSTMENTTIME~=0
    invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'adjustsampleclockrelativedelay',AWGADJUSTMENTTIME);
end

% Configure output mode - SEQ
invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_SEQ);
% % Configure output mode - ARB
% invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_ARB);

% Create arbitrary sequence

    % First part, line 2 to 511
    LineOfThePatternFirst = length(PatternInSampleFirst(:,1));
    SamplesOfOneLineFirst = length(PatternInSampleFirst(1,:));
    % Second part, line 512 to line 1
    SamplesOfOneLineLast = length(PatternInSampleLast);
    LineOfThePattern = LineOfThePatternFirst +1;
% WAVEFORMHANDLESARRAY = linspace(0,0,LineOfThePattern);
% First part, line 2 to 511
SAMPLECOUNTSARRAY=[];
for k=1:LineOfThePatternFirst
    WAVEFORMDATAARRAY = PatternInSampleFirst(k,:);
    WAVEFORMSIZE      = SamplesOfOneLineFirst;
    [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);
    WAVEFORMHANDLESARRAY(k)= NIFGEN_ATTR_ARB_WAVEFORM_HANDLE;
    SAMPLECOUNTSARRAY(k) = WAVEFORMSIZE;
end

% Second part, line 512 to line 1
    WAVEFORMDATAARRAY = PatternInSampleLast;
    WAVEFORMSIZE      = SamplesOfOneLineLast;
    [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);
    WAVEFORMHANDLESARRAY(LineOfThePattern)= NIFGEN_ATTR_ARB_WAVEFORM_HANDLE;
SAMPLECOUNTSARRAY(LineOfThePattern) = WAVEFORMSIZE;
    
 

SequenceLength = length(WAVEFORMHANDLESARRAY)
LOOPCOUNTSARRAY = linspace(1,1,SequenceLength);% Repeat time
% Advanced arb sequence
NIFGEN_VAL_NO_MARKER=-1;



MARKERLOCATIONARRAY=linspace(1,1,SequenceLength);
MARKERLOCATIONARRAY(1:SequenceLength)=NIFGEN_VAL_NO_MARKER;
MARKERLOCATIONARRAY(1) = 0;

COERCEDMARKERSARRAY=linspace(0,0,SequenceLength);
[COERCEDMARKERSARRAY,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'createadvancedarbsequence',SequenceLength,WAVEFORMHANDLESARRAY,LOOPCOUNTSARRAY,SAMPLECOUNTSARRAY,MARKERLOCATIONARRAY,COERCEDMARKERSARRAY);



% Query maximum capabilities
% [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities')
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);




% % Create arbitrary waveform
% [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);

% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);



% % Configure arbitrary waveform
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configurearbwaveform',ChannelName,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);
% Configure arbitrary sequency
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'configurearbsequence',ChannelName,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);



% Tigger mode :Continuouse/Stepped
if StepOrContinuous ==1
    invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuretriggermode',ChannelName,NIFGEN_VAL_STEPPED);
    invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuredigitaledgestarttrigger',"PFI0",Rising_Edge);
end


% export Marker 0 to PFI 1
invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'exportsignal',NIFGEN_VAL_MARKER_EVENT,"Marker0","PFI1");

% Continuous mode [The default mode is continuous mode]
% NIFGEN_VAL_OPERATE_CONTINUOUS=0;
% invoke(handles.AWGdeviceObj.Configuration,'configureoperationmode',ChannelName,NIFGEN_VAL_OPERATE_CONTINUOUS)

% Initiate the Waveform Generation
invoke(handles.AWGdeviceObj.Waveformcontrol,'initiategeneration');

% Enable the Output
invoke(handles.AWGdeviceObj.Configuration,'configureoutputenabled', ChannelName, Enabled);
% disp('AWG is working');
NewhObject=hObject;
Newhandles=handles;


% 3.1 Run last waveform - continue
function [NewhObject, Newhandles]=RunLastWaveform(hObject, handles,StepOrContinuous)
% The feedback loop function starts from here
fprintf('\n\n\n');
fprintf('Run last calculated waveform\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Status check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AWG Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice ~= 1
    errordlg('AWG is not intialized','GUI Error');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters define %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters define, these values are defined in the driver's header file 'niFgen.h'
NIFGEN_VAL_OUTPUT_ARB = 1;      % Arbitrary mode
NIFGEN_VAL_OUTPUT_SEQ = 2;
ChannelName = '0';              % This value is described in the help file 'NI Signal Generators Help'
Enabled = 1;                    % Enable output waveform
Rising_Edge = 101;              % Start trigger: rising edge detection
Falling_Edge = 102;             % Start trigger: falling edge detection
NIFGEN_VAL_CONTINUOUS = 2;      % Trigger mode: continuous 
NIFGEN_VAL_STEPPED = 3;         % Trigger mode: stepped 
NIFGEN_VAL_MARKER_EVENT = 1001; % Marker event, for Marker 0
% FirstTest =0;


%%% AWG group
% Sampling rate, E15, in Hz
AWGSamplingRate = str2double(get(handles.edit15, 'String'))*1e6;     
if AWGSamplingRate>114e6 || AWGSamplingRate<=0.01e6
    errordlg('Sampling range is [0.01, 114] MHz','GUI Error');
    return;
end
% Gain, E6, Range (0,6)
NIFGEN_ATTR_ARB_GAIN = str2double(get(handles.edit6, 'String'));
if NIFGEN_ATTR_ARB_GAIN>6 || NIFGEN_ATTR_ARB_GAIN<=0
    errordlg('Gain range is [0, 6]','GUI Error');
    return;
end
% Offset, E7, Range: (-0.25,0.25)
NIFGEN_ATTR_ARB_OFFSET = str2double(get(handles.edit7, 'String'));   
if NIFGEN_ATTR_ARB_OFFSET>0.25 || NIFGEN_ATTR_ARB_OFFSET<-0.25
    errordlg('AWG offset range is [-0.25, 0.25]','GUI Error');
    return;
end
% Sampling source, T26, P14, P15, 1. Internal 2. External
AWGSamplingSource = get(handles.text26, 'String');      
% AWG output delay, this value should be less than 1 sample time
AWGADJUSTMENTTIME = str2double(get(handles.edit43, 'String'))*1e-6; % Read how loops

WAVEFORMDATAARRAY = handles.LastWAVEFORMDATAARRAY;          % Waveform, value: (-1,1)
WAVEFORMSIZE = length(WAVEFORMDATAARRAY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Enable AWG output %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Disable the AWG if it is running. (To protect AWG)
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice == 1
    invoke(handles.AWGdeviceObj.Utility,'disable');
end



% [hObject,handles] = AWGInitializationFunction(hObject, 1, handles);
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'cleararbmemory');

% set(handles.AWGdeviceObj.Output(1), 'Output_Enabled', 0.0);
% invoke(handles.AWGdeviceObj.Waveformcontrol,'abortgeneration');
    

% if kLoopnum>1
% %     Configurationfunctionsarbitrarysequenceoutput
%     invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);
% end

% Sample rate multiplier: may be this is for output sample clock!!!
% set(handles.AWGdeviceObj.Clocksadvanced, 'External_Sample_Clock_Multiplier', 3)
% get83 = get(handles.AWGdeviceObj.Clocksadvanced, 'External_Sample_Clock_Multiplier')

% Configure the sampling rate
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configuresamplerate',AWGSamplingRate);
if AWGSamplingSource=='External' % The next line is only for external source.
   invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'configuresampleclocksource',"ClkIn");
end

% Configure AWG sample delay
if AWGADJUSTMENTTIME~=0
    invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'adjustsampleclockrelativedelay',AWGADJUSTMENTTIME);
end

% Configure output mode - SEQ
invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_SEQ);
% % Configure output mode - ARB
% invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_ARB);

% Create arbitrary sequence
[NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);
WAVEFORMHANDLESARRAY = [];

for k=1
    WAVEFORMHANDLESARRAY(k)=NIFGEN_ATTR_ARB_WAVEFORM_HANDLE;
end
SequenceLength = length(WAVEFORMHANDLESARRAY);
LOOPCOUNTSARRAY = linspace(1,1,SequenceLength);% Repeat time
% Advanced arb sequence
NIFGEN_VAL_NO_MARKER=-1;

SAMPLECOUNTSARRAY=linspace(1,1,SequenceLength);
SAMPLECOUNTSARRAY(1:SequenceLength) = WAVEFORMSIZE;
% SAMPLECOUNTSARRAY=[WAVEFORMSIZE,WAVEFORMSIZE];

MARKERLOCATIONARRAY=linspace(1,1,SequenceLength);
MARKERLOCATIONARRAY(1:SequenceLength)=NIFGEN_VAL_NO_MARKER;
MARKERLOCATIONARRAY(1) = 0;

COERCEDMARKERSARRAY=linspace(0,0,SequenceLength);
[COERCEDMARKERSARRAY,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'createadvancedarbsequence',SequenceLength,WAVEFORMHANDLESARRAY,LOOPCOUNTSARRAY,SAMPLECOUNTSARRAY,MARKERLOCATIONARRAY,COERCEDMARKERSARRAY);



% Query maximum capabilities
% [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities')
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);




% % Create arbitrary waveform
% [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);

% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'cleararbwaveform',NIFGEN_ATTR_ARB_WAVEFORM_HANDLE);



% % Configure arbitrary waveform
% invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configurearbwaveform',ChannelName,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);
% Configure arbitrary sequency
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'configurearbsequence',ChannelName,NIFGEN_ATTR_ARB_SEQUENCE_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);



% Tigger mode :Continuouse/Stepped
if StepOrContinuous ==1
    invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuretriggermode',ChannelName,NIFGEN_VAL_STEPPED);
    invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuredigitaledgestarttrigger',"PFI0",Rising_Edge);
end


% export Marker 0 to PFI 1
invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'exportsignal',NIFGEN_VAL_MARKER_EVENT,"Marker0","PFI1");

% Continuous mode [The default mode is continuous mode]
% NIFGEN_VAL_OPERATE_CONTINUOUS=0;
% invoke(handles.AWGdeviceObj.Configuration,'configureoperationmode',ChannelName,NIFGEN_VAL_OPERATE_CONTINUOUS)

% Initiate the Waveform Generation
invoke(handles.AWGdeviceObj.Waveformcontrol,'initiategeneration');

% Enable the Output
invoke(handles.AWGdeviceObj.Configuration,'configureoutputenabled', ChannelName, Enabled);
% disp('AWG is working');
NewhObject=hObject;
Newhandles=handles;

% 3.x Feedback loop main function - Run
function [NewhObject, Newhandles]=FeedbackloopRun(hObject, handles)
% The feedback loop function starts from here
fprintf('\n\n\n');
fprintf('Feedback loop starts from here.\n');
FeedbackLoopTime = tic;
FeedbackLoopTimeBeforeLoop = tic;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Status check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AWG Initialization status update
ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
if ExistAWGICdevice ~= 1
    errordlg('AWG is not intialized','GUI Error');
    return;
end


% Configure the board's sample rate, input, and trigger settings
if ~configureBoard(handles)
    fprintf('Error: DAQ/OSC Board configuration failed\n');
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters define %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters define, these values are defined in the driver's header file 'niFgen.h'
NIFGEN_VAL_OUTPUT_ARB = 1;      % Arbitrary mode
ChannelName = '0';              % This value is described in the help file 'NI Signal Generators Help'
Enabled = 1;                    % Enable output waveform
Rising_Edge = 101;              % Start trigger: rising edge detection
Falling_Edge = 102;             % Start trigger: falling edge detection
NIFGEN_VAL_CONTINUOUS = 2;      % Trigger mode: continuous 
NIFGEN_VAL_STEPPED = 3;         % Trigger mode: stepped 
NIFGEN_VAL_MARKER_EVENT = 1001; % Marker event, for Marker 0
% FirstTest =0;


%%% AWG group
% Sampling rate, E15, in Hz
AWGSamplingRate = str2double(get(handles.edit15, 'String'))*1e6;     
if AWGSamplingRate>114e6 || AWGSamplingRate<=0.01e6
    errordlg('Sampling range is [0.01, 114] MHz','GUI Error');
    return;
end
% Gain, E6, Range (0,6)
NIFGEN_ATTR_ARB_GAIN = str2double(get(handles.edit6, 'String'));
if NIFGEN_ATTR_ARB_GAIN>6 || NIFGEN_ATTR_ARB_GAIN<=0
    errordlg('Gain range is [0, 6]','GUI Error');
    return;
end
% Offset, E7, Range: (-0.25,0.25)
NIFGEN_ATTR_ARB_OFFSET = str2double(get(handles.edit7, 'String'));   
if NIFGEN_ATTR_ARB_OFFSET>0.25 || NIFGEN_ATTR_ARB_OFFSET<-0.25
    errordlg('AWG offset range is [-0.25, 0.25]','GUI Error');
    return;
end
% Sampling source, T26, P14, P15, 1. Internal 2. External
AWGSamplingSource = get(handles.text26, 'String');      
% AWG output delay, this value should be less than 1 sample time
AWGADJUSTMENTTIME = str2double(get(handles.edit43, 'String'))*1e-6; % Read how loops

%%% OSC group
% Acquiredata Start
OSCboardHandle = handles.boardHandle;
% call mfile with library definitions
AlazarDefs;
% Sampling rate set, edit 11, in Hz
OSCSamplingRate = str2double(get(handles.edit11, 'String'))*1e6;
if OSCSamplingRate~=500e6 && OSCSamplingRate~=1000e6 && OSCSamplingRate~= 912e6
    errordlg('OSC sampling rate should be 500 or 1000 or 912 MHz','GUI Error');
    return;
end
% Buffer time, edit 16, in second
OSCBufferTime = str2double(get(handles.edit16, 'String'))*1e-6;
% Buffer samples
OSCBufferSamples = OSCBufferTime*OSCSamplingRate;
% Buffer numbers in 1 acquisition
OSCBufferNumberInOneAcq = str2double(get(handles.edit12, 'String'));
% Acquisition Time, in second
OSCAcquisitionTime = OSCBufferTime * OSCBufferNumberInOneAcq;
% Samples in one total acquisition
OSCNlengthAllChannel = OSCAcquisitionTime*OSCSamplingRate;
% Channels
OSCChannelNumber = 1;
OSCChannelNumberCode = get(handles.text50, 'String');
if OSCChannelNumberCode == 'A+B'
    OSCChannelNumber = 2;
end



% Samples in one channel in one acquisition, last 3 samples are deleted.
OSCNlengthPerChannel = floor(OSCNlengthAllChannel)-3;
% OSC time
DTimeOSC = 1/OSCSamplingRate;

VTimeOSC = 1:1:OSCNlengthPerChannel;
VTimeOSC = VTimeOSC.*DTimeOSC;
% VTimeOSC = linspace(1,OSCNlengthPerChannel,OSCNlengthPerChannel);
% VTimeOSC = DTimeOSC* VTimeOSC;

% Ratio of OSC and AWG
SamRatioOscToAWG = round(OSCSamplingRate/AWGSamplingRate);              % Ratio of OSC sampling rate and AWG sampling rate


%%% Feedback loop group
NLoop = str2double(get(handles.edit17, 'String')); % Read how loops
DataSaveWorkModeNum = 666;
% DataSaveWorkMode = get(handles.text24, 'String'); % 1. Test mode 2. Work mode
% Show Raw data or not?
FeedbackLoopRawData = get(handles.text108, 'String'); % 1. Test mode 2. Work mode

% Data Save Group
DataSaveSaveData = get(handles.text25, 'String'); % 1. Not save data 2. Save data

% Save data - check if same name exist
if DataSaveSaveData == 'Yes'
    PresentFolder = get(handles.edit23, 'String');
    PresentFileName = [get(handles.edit1, 'String'),'_',get(handles.edit2, 'String'),'_',get(handles.edit3, 'String')];
    PresentAllFilePath = [PresentFolder,'\',PresentFileName,'.mat'];
    y=exist(PresentAllFilePath,'file');
    if y == 2
    message = sprintf('The file already exsited, cover it or not?');
	reply = questdlg(message, 'Cover present file', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		return;
    end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% DAQ/OSC predefine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Select which channels to capture (A, B, or both)

% TODO: Select which channels to capture (A, B, or both)
OSCChannel = get(handles.text50, 'String');
if OSCChannel == 'A  '
    OSCchannelMask = CHANNEL_A;
end
if OSCChannel == 'B  '
    OSCchannelMask = CHANNEL_B;
    SamplingClockSource = get(handles.text134, 'String'); % Internal or External
    if SamplingClockSource == 'External'
        fprintf('Only CH A is available when use external source as sampling clock.\n', errorToText(retCode));
        return
    end
end
if OSCChannel == 'A+B'
    OSCchannelMask = CHANNEL_A+CHANNEL_B;
    SamplingClockSource = get(handles.text134, 'String'); % Internal or External
%     % 20220606 test
%     if SamplingClockSource == 'External'
%         fprintf('Only CH A is available when use external source as sampling clock.\n', errorToText(retCode));
%         return
%     end
end

% Calculate the number of enabled channels from the channel mask
OSCchannelCount = 0;
OSCchannelsPerBoard = 2;
for channel = 0:OSCchannelsPerBoard - 1
    channelId = 2^channel;
    if bitand(channelId, OSCchannelMask)
        OSCchannelCount = OSCchannelCount + 1;
    end
end

if (OSCchannelCount < 1) || (OSCchannelCount > OSCchannelsPerBoard)
    fprintf('Error: Invalid channel mask %08X\n', OSCchannelMask);
    return
end
% Get the sample and memory size
[OSCretCode, OSCboardHandle, maxSamplesPerRecord, OSCbitsPerSample] = AlazarGetChannelInfo(OSCboardHandle, 0, 0);
if OSCretCode ~= ApiSuccess
    fprintf('Error: AlazarGetChannelInfo failed -- %s\n', errorToText(OSCretCode));
    return
end
% Calculate the size of each buffer in bytes
OSCbytesPerSample = floor((double(OSCbitsPerSample) + 7) / double(8));
OSCsamplesPerBuffer = OSCBufferSamples * OSCchannelCount;
OSCbytesPerBuffer = OSCbytesPerSample * OSCsamplesPerBuffer;
% Find the number of buffers in the acquisition
if OSCAcquisitionTime > 0
    OSCsamplesPerAcquisition = floor((OSCSamplingRate * OSCAcquisitionTime + 0.5));
    OSCbuffersPerAcquisition = uint32(floor((OSCsamplesPerAcquisition + OSCBufferSamples - 1) / OSCBufferSamples));
else
    OSCbuffersPerAcquisition = hex2dec('7FFFFFFF');  % acquire until aborted
end
% TODO: Select the number of DMA buffers to allocate.
% The number of DMA buffers must be greater than 2 to allow a board to DMA into
% one buffer while, at the same time, your application processes another buffer.
OSCbufferCount = uint32(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% AWG waveform define %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PulseOn = 1;                % Waveform value when pulse is on 
PulseOff = 0;               % Waveform value when pulse is off
NIFGEN_VAL_OPERATE_CONTINUOUS=0; 
DTimeAWG = 1/AWGSamplingRate;      % s
Scan_Comp_Sample      = str2double(get(handles.edit51, 'String'));
ScanBackSamples_left     = str2double(get(handles.edit47, 'String')); %Scan back time for left
ScanSamples              = str2double(get(handles.edit48, 'String'));
ScanBackSamples_Right    = str2double(get(handles.edit49, 'String'));
ScanBackSamples_Right2   = str2double(get(handles.edit50, 'String'));
FlyBackSamples        = str2double(get(handles.edit52, 'String')); %Flyback time, in samples of AWG.
% NSampleAWG            = floor(ScanAllLineSamples/4)*4; % 3772 
% VIntensityAWG         = linspace(PulseOff,PulseOff,NSampleAWG); % Pre-dine AWG pattern
% VIntensityAWG(560+1)  = PulseOn;
% VIntensityAWG(600:1000) = PulseOn;
% VIntensityAWG(2600:2700) = PulseOn;
% WAVEFORMDATAARRAY     = VIntensityAWG;
% VTimeAWG              = linspace(DTimeAWG,DTimeAWG*NSampleAWG,NSampleAWG);          % Define Time axis


    
% if image is converted to pattern
if get(handles.text94, 'String') == 'Yes'

    % First part, line 2 to 511
    PatternInSampleFirst = handles.ConvertedPatternInSamplesFirst;
    LineOfThePatternFirst = length(PatternInSampleFirst(:,1)); % number of lines
    SamplesOfOneLineFirst = length(PatternInSampleFirst(1,:)); % samples of a line
    SamplesOfOnePatternFirst = (SamplesOfOneLineFirst+Scan_Comp_Sample)*LineOfThePatternFirst; % total samples
    % Second part, line 512 to line 1
    PatternInSampleLast = handles.ConvertedPatternInSamplesLast;
    SamplesOfOneLineLast = length(PatternInSampleLast);
    
    SamplesOfOnePattern = SamplesOfOnePatternFirst + SamplesOfOneLineLast;
    SamplesOfOnePattern=floor(SamplesOfOnePattern/4)*4+4;
    NSampleAWG = SamplesOfOnePattern;
    WAVEFORMDATAARRAY = linspace(0,0,SamplesOfOnePattern);
    WAVEFORMDATAARRAYImageArea = WAVEFORMDATAARRAY;
    
    
    % First part, line 2 to 511
    for k = 1:LineOfThePatternFirst
        StartL = 1+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
        EndL   = SamplesOfOneLineFirst+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
        WAVEFORMDATAARRAY(StartL:EndL)=PatternInSampleFirst(k,:);
        WAVEFORMDATAARRAYImageArea(StartL:EndL)=handles.ConvertedImageArea(k,:);
    end
    
    % Second part, line 512 to 1
    StartL = 1+SamplesOfOnePatternFirst;
    EndL   = SamplesOfOneLineLast+SamplesOfOnePatternFirst;
    WAVEFORMDATAARRAY(StartL:EndL)=PatternInSampleLast(:);
    
%     NSampleAWG = 2*NSampleAWG;
%     WAVEFORMDATAARRAY = [WAVEFORMDATAARRAY,WAVEFORMDATAARRAY];
    WAVEFORMDATAARRAY = WAVEFORMDATAARRAY/1*(PulseOn-PulseOff)+PulseOff; % from [0,1] to [-1,1]    
    VTimeAWG = linspace(DTimeAWG,DTimeAWG*NSampleAWG,NSampleAWG);          % Define Time axis
  
    
end

    FigXPeriod       = (str2double(get(handles.edit57, 'String'))-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample)*DTimeAWG; % every two lines have 2606 AWG samples.
    FigXStartValue = str2double(get(handles.edit58, 'String'))*1e-6+FigXPeriod;
    FigXEndValue   = str2double(get(handles.edit56, 'String'))*1e-6+FigXPeriod;
    FigureShowStart = find(VTimeAWG>FigXStartValue,1);
    FigureShowEnd = find(VTimeAWG>FigXEndValue,1);

    
    handles.LastWAVEFORMDATAARRAY     = WAVEFORMDATAARRAY;
    handles.FinalPatternInSampleFirst = PatternInSampleFirst;
    handles.FinalPatternInSampleLast = PatternInSampleLast;
    
% export initial AWG waveform
csvwrite('patternB4Fdbk.csv', WAVEFORMDATAARRAY);

% WAVEFORMDATAARRAY = WAVEFORMDATAARRAY.*0+1;
% Line = 772
% How many period does one acquisition has?*
% AverageTimeCalc = floor(OSCAcquisitionTime/TWindowAWG);
% set(handles.edit28, 'String',num2str(AverageTimeCalc));

WAVEFORMDATAARRAYNorm = (WAVEFORMDATAARRAY-PulseOff)/(PulseOn-PulseOff);
WAVEFORMSIZE = length(WAVEFORMDATAARRAY);               % Waveform size
% WaveformIdealOut = WAVEFORMDATAARRAY;                   % Ideal final optical output waveform, normalized.
WaveformIdealMatrix = reshape(WAVEFORMDATAARRAY,[1,WAVEFORMSIZE/1]); % Matrix
WaveformIdealPeak = (WaveformIdealMatrix-PulseOff)/(PulseOn-PulseOff);% Ideal final output peak.
WavefommIdealPulseOnOne = WaveformIdealPeak>=PulseOn;   % Pulse that is >=1
WavefommIdealPulseOn    = WaveformIdealPeak>PulseOff;          % Pulses that is >0
NumPulseOnOne = sum(WavefommIdealPulseOnOne);                    % number of pulse on
NumPulseOn    = sum(WavefommIdealPulseOn);                    % number of pulse on

% figure
%
%plot(VTimeAWG,WAVEFORMDATAARRAY,'b');hold on
%plot(VTimeAWG,WaveformIdealPeak,'r');hold on
%plot(VTimeAWG,WavefommIdealPulseOnOne,'--k');hold on
%plot(VTimeAWG,WavefommIdealPulseOn,'--g');hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Feedback loop main function %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% record video? 1 for yes, 0 for no. 
videoRecord = 0; 

if ishandle(211)
    close figure 211;
end
if ishandle(212)
    close figure 212;
end
if ishandle(213)
    close figure 213;
end

fh = figure (211);
fh.Units = 'normalized'; fh.Position = [0 0 1 1];          % maximize figure 11 when run 
hold on
% save(['Feedback_loop',Date,'_',ObjectName,'_',GroupName,'_',num2str(k),'th_rounds'],'VTime','VIntAmpPeakArrayOn','VTimeOSC','VItensityAmp');
if FeedbackLoopRawData == 'Yes'
    set(gcf, 'units','normalized','outerposition',[0 0.035 0.5 0.7]);
    figure (212)
    set(gcf, 'units','normalized','outerposition',[0.5 0.035 0.5 0.7]);
else
    set(gcf, 'units','normalized','outerposition',[0 0.035 0.5 0.7]);
end



% Arrays define
VIntAmpPeakArrayMeanOn = linspace(0,0,NLoop);             % Average output pulse peak
VIntAmpPeakArrayOnAveFluc =  VIntAmpPeakArrayMeanOn;       % average fluctuation
VIntAmpPeakArrayOnRMSFluc =  VIntAmpPeakArrayMeanOn;       % RMS fluctuation
VIntAmpPeakArrayOnMaxFluc =  VIntAmpPeakArrayMeanOn;       % Max fulctuatopm
VIntAmpPeakArrayOnOneAveFluc =  VIntAmpPeakArrayMeanOn;       % average fluctuation
VIntAmpPeakArrayOnOneRMSFluc =  VIntAmpPeakArrayMeanOn;       % RMS fluctuation
VIntAmpPeakArrayOnOneMSFluc =  VIntAmpPeakArrayMeanOn;       % Mean squared fluctuation
VIntAmpPeakArrayOnOneMaxFluc =  VIntAmpPeakArrayMeanOn;       % Max fulctuatopm
VLoops = linspace(1,NLoop,NLoop);                       % Array of loops






set(handles.AWGdeviceObj.Output(1), 'Output_Enabled', 0.0);
invoke(handles.AWGdeviceObj.Waveformcontrol,'abortgeneration');
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'cleararbmemory');

% Configure the sampling rate
invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configuresamplerate',AWGSamplingRate);
if AWGSamplingSource=='External' % The next line is only for external source.
   invoke(handles.AWGdeviceObj.Configurationfunctionsconfigureclock,'configuresampleclocksource',"ClkIn");
end

PulsePeakTest = str2double(get(handles.edit53, 'String'));



    % Configure the board's sample rate, input, and trigger settings
if ~configureBoard(handles)
    fprintf('Error: Board configuration failed\n');
    return
end

% Feedback loop: Acquire data, optionally saving it to a file
AverageTime = str2double(get(handles.edit28, 'String'));
VItensityAmp = VTimeOSC*0;
VItensityMod = VItensityAmp;

% Set the size of the waveform (equal to AWG waveform after scaling)
WaveformOSCSize = round(WAVEFORMSIZE*SamRatioOscToAWG);                
% Find when does the waveform start
ChannelAStartPoint = str2double(get(handles.edit29, 'String'))*1e-6;
WaveChannelAStartPoint = find(VTimeOSC>=ChannelAStartPoint & VTimeOSC<(ChannelAStartPoint)+0.01e-6); % Find when does the waveform start
if WaveChannelAStartPoint(1)+WaveformOSCSize-1>OSCNlengthPerChannel
    errordlg('OSC record length is less than the AWG waveform period','GUI Error');
    return;
end

ChannelBStartPoint = str2double(get(handles.edit30, 'String'))*1e-6;
WaveChannelBStartPoint = find(VTimeOSC>ChannelBStartPoint & VTimeOSC<(ChannelBStartPoint)+0.01e-6); % Find when does the waveform start
    

% Amplified intensity signal, the start and end
VItensityAmpPulseStart = WaveChannelAStartPoint(1);
VItensityAmpPulseEnd   = VItensityAmpPulseStart + SamRatioOscToAWG * WAVEFORMSIZE -1;
disp(int2str(VItensityAmpPulseStart))

% Modulated intensity signal, the start and end
VItensityModPulseStart = WaveChannelBStartPoint(1);
VItensityModPulseEnd   = VItensityModPulseStart + SamRatioOscToAWG * WAVEFORMSIZE -1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWG pre-setting %%%%%%%%%%%%%%%%%%%%%%%%
% export Marker 0 to PFI 1, 0.04 s
set(handles.AWGdeviceObj.Arbitrarywaveformarbitrarywaveformmode(1), 'Marker_Position', 0);
invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'exportsignal',NIFGEN_VAL_MARKER_EVENT,"Marker0","PFI1");

% Configure output mode, 0.01 s
invoke(handles.AWGdeviceObj.Configuration,'configureoutputmode',NIFGEN_VAL_OUTPUT_ARB);

% Tigger mode :Continuouse/Stepped
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuretriggermode',ChannelName,NIFGEN_VAL_STEPPED);
% invoke(handles.AWGdeviceObj.Configurationfunctionstriggeringandsynchronization,'configuredigitaledgestarttrigger',"PFI0",Rising_Edge);

% Continuous mode [The default mode is continuous mode], 0.01 s
% invoke(handles.AWGdeviceObj.Configuration,'configureoperationmode',ChannelName,NIFGEN_VAL_OPERATE_CONTINUOUS);
PulsePeakTest = str2double(get(handles.edit53, 'String'));
fprintf(['Before loop: ',num2str(toc(FeedbackLoopTimeBeforeLoop)),' s.\n']);

%csvwrite('waveformIdealPulseOnOne.csv', WavefommIdealPulseOnOne);
%csvwrite('waveformDataArray.csv', WAVEFORMDATAARRAY);
% csvwrite('waveformDataArrayNorm.csv', WAVEFORMDATAARRAYNorm);
% csvwrite('waveformIdealPulseOn.csv', WavefommIdealPulseOn);
% csvwrite('waveformSize.csv', WAVEFORMSIZE);

lastRMS = 100;                          % used to track the most recent RMS, so skip the current loop if RMS jumps by too much 

% Feedback loop starts
for kLoopnum=1:NLoop
    WAVEFORMDATAARRAY(WAVEFORMSIZE) = 0; % for safty reason, so that the modulated light alway has some power. It takes 0.0023903 seconds
    
    AWGTakesTime = tic;
    % WAVEFORMSIZE = length(WAVEFORMDATAARRAY); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Enable AWG output %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disable the AWG if it is running. (To protect AWG)
    % ExistAWGICdevice = isfield(handles,'AWGdeviceObj');
    % if ExistAWGICdevice == 1
    %     invoke(handles.AWGdeviceObj.Utility,'disable');
    % end
    % [hObject,handles] = AWGInitializationFunction(hObject, 1, handles);

%     whos WAVEFORMDATAARRAY

    % clear AWG memory so it won't be loaded with unneeded waveforms from
    % previous loops -- by SZhao 01062021
    invoke(handles.AWGdeviceObj.Waveformcontrol,'abortgeneration');
    invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarysequenceoutput,'cleararbmemory');
        
    % Create arbitrary waveform, 0.074 s .... 
%                                230 ms for 3.8 M pulses
    [NIFGEN_ATTR_ARB_WAVEFORM_HANDLE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'createwaveformf64',ChannelName,WAVEFORMSIZE,WAVEFORMDATAARRAY);


    % Configure arbitrary waveform, 0.01 s
    invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'configurearbwaveform',ChannelName,NIFGEN_ATTR_ARB_WAVEFORM_HANDLE,NIFGEN_ATTR_ARB_GAIN,NIFGEN_ATTR_ARB_OFFSET);

    % Initiate the Waveform Generation, 0.016 s
    invoke(handles.AWGdeviceObj.Waveformcontrol,'initiategeneration');

    % Enable the Output, 0.02 s
    invoke(handles.AWGdeviceObj.Configuration,'configureoutputenabled', ChannelName, Enabled);
    % disp('AWG is working');
    fprintf(['AWG: ',num2str(round(toc(AWGTakesTime)*1000)),' ms.   ']);
    
    pause(0.2);               % by SZhao 20210823, to let the DAQ after PD wait so that the transient from IM has settled down
    
    % Query maximum capabilities
%     [MAXIMUMNUMBEROFWAVEFORMS,WAVEFORMQUANTUM,MINIMUMWAVEFORMSIZE,MAXIMUMWAVEFORMSIZE] = invoke(handles.AWGdeviceObj.Configurationfunctionsarbitrarywaveformoutput,'queryarbwfmcapabilities');
%     disp('Max number of waveforms is ' + string(MAXIMUMNUMBEROFWAVEFORMS) + 'Max waveform size is ' + string(MAXIMUMWAVEFORMSIZE));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Read OSC waveform %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acquiredata starts
    OSCTakesTime = tic;
    if AverageTime~=1 % if average time is not 1, then these arrays should be reset.
        VItensityAmp = VItensityAmp*0;
        if OSCChannelNumber == 2
            VItensityMod = VItensityMod*0;
        end
    end

    % Acquire several data to do average
    for AverageNum=1:AverageTime
        % Configure the board's sample rate, input, and trigger settings
        % if ~configureBoard(handles)
        %     fprintf('Error: Board configuration failed\n');
        %     return
        % end
        
        % set default return code to indicate failure
        result = false;
        % Create an array of DMA buffers
        OSCbuffers = cell(1, OSCbufferCount);
        for j = 1 : OSCbufferCount
            OSCpbuffer = AlazarAllocBuffer(OSCboardHandle, OSCbytesPerBuffer);
            if OSCpbuffer == 0
                fprintf('Error: AlazarAllocBuffer %u samples failed\n', OSCsamplesPerBuffer);
                return
            end
            OSCbuffers(1, j) = { OSCpbuffer };
        end

        OSCretCode = AlazarBeforeAsyncRead(OSCboardHandle, OSCchannelMask, 0, OSCBufferSamples, 1, hex2dec('7FFFFFFF'), ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRIGGERED_STREAMING);
        if OSCretCode ~= ApiSuccess
            fprintf('Error: AlazarBeforeAsyncRead failed -- %s\n', errorToText(OSCretCode));
            return
        end
        % Post the buffers to the board
        for bufferIndex = 1 : OSCbufferCount
            OSCpbuffer = OSCbuffers{1, bufferIndex};
            OSCretCode = AlazarPostAsyncBuffer(OSCboardHandle, OSCpbuffer, OSCbytesPerBuffer);
            if OSCretCode ~= ApiSuccess
                fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(OSCretCode));
                return
            end
        end
        % Arm the board system to wait for triggers
        OSCretCode = AlazarStartCapture(OSCboardHandle);
        if OSCretCode ~= ApiSuccess
            fprintf('Error: AlazarStartCapture failed -- %s\n', errorToText(OSCretCode));
            return
        end

        OSCbuffersCompleted = 0;
        OSCcaptureDone = false;
        OSCsuccess = false;
        % if more than one buffer is used, we need to pre-define OSCFinaldata.
        if OSCBufferNumberInOneAcq ~=1 
            OSCFinaldata =[];
        end


        % OSC working ... 
        while ~OSCcaptureDone
            OSCbufferIndex = mod(OSCbuffersCompleted, OSCbufferCount) + 1;
            OSCpbuffer = OSCbuffers{1, OSCbufferIndex};

            % Wait for the first available buffer to be filled by the board
            [OSCretCode, OSCboardHandle, OSCbufferOut] = ...
                AlazarWaitAsyncBufferComplete(OSCboardHandle, OSCpbuffer, 5000);
            if OSCretCode == ApiSuccess
                % This buffer is full
                OSCbufferFull = true;
                OSCcaptureDone = false;
            elseif OSCretCode == ApiWaitTimeout
                % The wait timeout expired before this buffer was filled.
                % The board may not be triggering, or the timeout period may be too short.
                fprintf('Error: AlazarWaitAsyncBufferComplete timeout -- Verify trigger!\n');
                OSCbufferFull = false;
                OSCcaptureDone = true;
            else
                % The acquisition failed
                fprintf('Error: AlazarWaitAsyncBufferComplete failed -- %s\n', errorToText(OSCretCode));
                OSCbufferFull = false;
                OSCcaptureDone = true;
            end

            if OSCbufferFull
                % TODO: Process sample data in this buffer.
                %
                % NOTE:
                %
                % While you are processing this buffer, the board is already
                % filling the next available buffer(s).
                %
                % You MUST finish processing this buffer and post it back to the
                % board before the board fills all of its available DMA buffers
                % and on-board memory.
                %
                % Records are arranged in the buffer as follows: R0A, R1A, R2A ... RnA, R0B,
                % R1B, R2B ...
                % with RXY the record number X of channel Y
                %
                % A 12-bit sample code is stored in the most significant bits of
                % in each 16-bit sample value.
                %
                % Sample codes are unsigned by default. As a result:
                % - a sample code of 0x0000 represents a negative full scale input signal.
                % - a sample code of 0x8000 represents a ~0V signal.
                % - a sample code of 0xFFFF represents a positive full scale input signal.
                if OSCbytesPerSample == 1
                    setdatatype(OSCbufferOut, 'uint8Ptr', 1, OSCsamplesPerBuffer);
                else
                    setdatatype(OSCbufferOut, 'uint16Ptr', 1, OSCsamplesPerBuffer);
                end
                % Save data
                if OSCBufferNumberInOneAcq ==1
                    OSCFinaldata = double(OSCbufferOut.Value);
                else
                    OSCFinaldata = [OSCFinaldata,double(OSCbufferOut.Value)];
                end
                % Make the buffer available to be filled again by the board
                OSCretCode = AlazarPostAsyncBuffer(OSCboardHandle, OSCpbuffer, OSCbytesPerBuffer);
                if OSCretCode ~= ApiSuccess
                    fprintf('Error: AlazarPostAsyncBuffer failed -- %s\n', errorToText(OSCretCode));
                    OSCcaptureDone = true;
                end
                % Update progress
                OSCbuffersCompleted = OSCbuffersCompleted + 1;
                if OSCbuffersCompleted >= OSCbuffersPerAcquisition
                    OSCcaptureDone = true;
                    OSCsuccess = true;
                end
            end % if OSCbufferFull
        end % while ~OSCcaptureDone

        % Abort the acquisition
        OSCretCode = AlazarAbortAsyncRead(OSCboardHandle);
        if OSCretCode ~= ApiSuccess
            fprintf('Error: AlazarAbortAsyncRead failed -- %s\n', errorToText(OSCretCode));
        end
        % Release the buffers
        for OSCbufferIndex = 1:OSCbufferCount
            OSCpbuffer = OSCbuffers{1, OSCbufferIndex};
            OSCretCode = AlazarFreeBuffer(OSCboardHandle, OSCpbuffer);
            if OSCretCode ~= ApiSuccess
                fprintf('Error: AlazarFreeBuffer failed -- %s\n', errorToText(OSCretCode));
            end
            clear OSCpbuffer;
        end
        % set return code to indicate OSCsuccess
        result = OSCsuccess;
        OSCFinaldata = (OSCFinaldata - 32768)*400/32768; % convert 2^16 data to mV.

        if result ~=1
            fprintf('Error: Acquisition failed\n');
            return
        end

        if OSCChannelNumber == 1
            if AverageTime==1
                VItensityAmp = OSCFinaldata(1:OSCNlengthPerChannel); % Channel A, amplified signal
            else
                VItensityAmp = VItensityAmp+OSCFinaldata(1:OSCNlengthPerChannel); % Channel A, amplified signal
            end
        else
            VItensityAmp = OSCFinaldata(1:2:(OSCNlengthPerChannel-1)*2+1); % Channel A, amplified signal
            VItensityMod = OSCFinaldata(2:2:(OSCNlengthPerChannel)*2); % Channel B, modulated signal
            
% %             test 20220607
%             VItensityAmp = VItensityAmp+OSCFinaldata(1:2:(OSCNlengthPerChannel-1)*2+1); % Channel A, amplified signal
%             VItensityMod = VItensityMod+OSCFinaldata(2:2:(OSCNlengthPerChannel)*2); % Channel B, modulated signal
        end
    end % End for average of OSC acquisition
    % Acquiredata End here

    fprintf(['OSC: ',num2str(round(toc(OSCTakesTime)*1000)),' ms.   ']);
    % figure(100)
    % plot(VItensityAmp);
    
%     invoke(handles.AWGdeviceObj.Waveformcontrol,'abortgeneration');  % for troubleshooting jumpiness, by SZhao 20210127
    PostProcessTakesTime = tic;

    % Depends on how many channels, there are two codes here. One channel code 
    % is in use and it is ready. Two channel code is not in use, so it is not optimized.
%     if OSCChannelNumber == 1       % used to match the end in line2896
        if AverageTime~=1 % average
            VItensityAmp = VItensityAmp/AverageTime;
            VItensityMod = VItensityMod/AverageTime;
        end
        % Save data
        if DataSaveSaveData == 'Yes'
            save([get(handles.edit23, 'String'),'\',PresentFileName,'_Amp_signal_Loop_',num2str(kLoopnum)],'VItensityAmp');
        end
        if FeedbackLoopRawData == 'Yes'
            figure (12)
            plot(VTimeOSC/1e-6,VItensityAmp);
            % hold on;
            xlabel('Time/us');
            ylabel('Voltage/mV');
            xlim([0,max(VTimeOSC/1e-6)-0.0001]);
            ylim([-100,400]);
            title('Channel A, amplified signal');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% processing waveform %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here assuming AWG sampling rate 32 MHz, laser rep rate 32 MHz, OSC sampling rate 672 MHz.
%         testTime = tic;

%         VIntAmpMatrix = vec2mat(VItensityAmp(VItensityAmpPulseStart:VItensityAmpPulseEnd),SamRatioOscToAWG); % 0.1956 s ... 0.26 s for 114 MHz laser    ******This is slow! Use reshape() instead.      
%         VIntAmpPeakArray = max(VIntAmpMatrix, [], 2); % Find the peak of each pulse
%         VIntAmpPeakArray = VIntAmpPeakArray'; % Take the transpose, so now VIntAmpPeakArray is a row vector (1 * n)
        
        % Use reshape() function to isolate each pulse into a column
        VIntAmpMatrix = reshape(VItensityAmp(VItensityAmpPulseStart:VItensityAmpPulseEnd),SamRatioOscToAWG,[]); % 8 * col matrix
        VIntAmpPeakArray = max(VIntAmpMatrix); % Find the peak of each pulse (find the max of each column) % output is 1 * n
        
        VIntModMatrix = reshape(VItensityMod(VItensityModPulseStart:VItensityModPulseEnd),SamRatioOscToAWG,[]); % 8 * col matrix
        VIntModPeakArray = max(VIntModMatrix); % Find the peak of each pulse (find the max of each column) % output is 1 * n
%         fprintf(['this step takes ',num2str(round(toc(testTime)*1000)),' ms.   \n'])  

%         whos VIntAmpMatrix
%         whos VIntAmpPeakArray
        if PulsePeakTest ~=0
            figure
            PulseStartTest = WaveChannelAStartPoint(1) + SamRatioOscToAWG * (PulsePeakTest-1);
            PulseEndTest   = PulseStartTest + SamRatioOscToAWG * 1 -1;
            plot(VTimeOSC(PulseStartTest:PulseEndTest)/1e-6,VItensityAmp(PulseStartTest:PulseEndTest),'o');
            hold on
            PulseStartTestWave = WaveChannelAStartPoint(1) + SamRatioOscToAWG * (PulsePeakTest-1-3);
            PulseEndTestWave   = PulseStartTestWave + SamRatioOscToAWG * 7 -1;
            plot(VTimeOSC(PulseStartTestWave:PulseEndTestWave)/1e-6,VItensityAmp(PulseStartTestWave:PulseEndTestWave),'-');
            hold on
            title(['Pulse ',num2str(PulsePeakTest)]);
            xlabel('Time/us');
            ylabel('Intensity/mV');
        end

        % Process pulses that are One!!! It takes 0.0161 seconds
        VIntAmpPeakArrayOnOne = VIntAmpPeakArray.*WavefommIdealPulseOnOne;           % Choose only area with pulse on
        VIntAmpPeakArrayOneOneMean(kLoopnum) = sum(VIntAmpPeakArrayOnOne)/NumPulseOnOne;   % Average value
        VIntAmpPeakArrayOnOneNorm = VIntAmpPeakArrayOnOne/max(VIntAmpPeakArrayOnOne);     % Normalize, [0, 1]
        VIntAmpPeakArrayOnOneNormMean = sum(VIntAmpPeakArrayOnOneNorm)/NumPulseOnOne;       % Mean value of pulse peak
        VIntAmpPeakArrayOnOneFluc = abs(VIntAmpPeakArrayOnOneNorm-VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
        VIntAmpPeakArrayOnOneFlucRatio = VIntAmpPeakArrayOnOneFluc/VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne;   % Fluctuation ratio
        VIntAmpPeakArrayOnOneAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio)/NumPulseOnOne;            % average fluctuation
        VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) = sqrt(sum(VIntAmpPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne);         % RMS fluctuation modified by SZhao 20210106
        VIntAmpPeakArrayOnOneMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne;         % Mean squared fluctuation modified by SZhao 20210106
        VIntAmpPeakArrayOnOneMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnOneFlucRatio);
        
%         % code below by SZhao 20210205 to determine noise floor
%         VIntAmpPeakArrayOnOneFluc2 = (VIntAmpPeakArrayOnOneNorm-VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
%         VIntAmpPeakArrayOnOneFlucRatio2 = VIntAmpPeakArrayOnOneFluc2/VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne;   % Fluctuation ratio
%         [maxRatio,indexMax] = max(VIntAmpPeakArrayOnOneFlucRatio2);
%         [minRatio,indexMin] = min(VIntAmpPeakArrayOnOneFlucRatio2);
%         disp('indices of max/min are ' + string(indexMax) + '  ' + string(indexMin));

        % Process pulses that are not off, It takes 0.0161 seconds 
        %                                           50 ms for 3.8 M pulses
        % can be discarded maybe to save 50 ms? 20210519
        VIntAmpPeakArrayOn = VIntAmpPeakArray.*WavefommIdealPulseOn;           % Choose only area with pulse on
        VIntAmpPeakArrayMeanOn(kLoopnum) = sum(VIntAmpPeakArrayOn)/NumPulseOn;   % Average value
        VIntAmpPeakArrayOnNorm = VIntAmpPeakArrayOn/max(VIntAmpPeakArrayOn);     % Normalize, [0, 1]
        VIntAmpPeakArrayOnNormMean = sum(VIntAmpPeakArrayOnNorm)/NumPulseOn;       % Mean value of pulse peak
        VIntAmpPeakArrayOnFluc = abs(VIntAmpPeakArrayOnNorm-VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak); % Deviation from ideal peak

        VIntAmpPeakArrayOnFlucRatio = VIntAmpPeakArrayOnFluc/VIntAmpPeakArrayOnNormMean;   % Fluctuation ratio
        VIntAmpPeakArrayOnAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio)/NumPulseOn;            % average fluctuation
        VIntAmpPeakArrayOnRMSFluc(kLoopnum) = sqrt(sum(VIntAmpPeakArrayOnFlucRatio.^2)/NumPulseOn);         % RMS fluctuation modified by SZhao 01062021
        VIntAmpPeakArrayOnMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnFlucRatio);
        
        if VIntAmpPeakArrayOnRMSFluc(kLoopnum) > 1.5*lastRMS                            % if the RMS is jumpy due to AWG unstable DC level, then skip this iteration
            continue;
        end
        lastRMS = VIntAmpPeakArrayOnRMSFluc(kLoopnum);

        % figure
        % plot(VTimeAWG,VIntAmpPeakArrayOnNorm,'b');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak,'k');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnFluc,'r');hold on
        % plot(VTimeAWG,WaveformIdealPeak,'--g');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne,'k');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnOneFluc,'r');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOn,'--k');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnOneNorm,'--g');hold on

        % figure
        % plot(VTimeOSC/1e-6-ChannelAStartPoint/1e-6,VItensityAmp);hold on;
        % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm*20,'k');hold on

        % hold on;
        % xlabel('Time/us');
        % ylabel('Voltage/mV');
        % 
        % fprintf(['Test 2 takes',num2str(toc(FeedbackLoopTimeBeforeLoop)),' seconds.\n']);
        % % 
        % figure (13)
        % plot(VTimeAWG/1e-6,VTimeAWG.*0+VIntAmpPeakArrayOnNormMean,'--r');hold on
        % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnNorm,'-.k');hold on
        % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnFluc,'-b');hold on
        % xlim([0,40]);
        
        

    %     subplot(2, 3, 4);
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnNorm);hold on
    %     if kLoopnum == NLoop
    %         % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %         plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*VIntAmpPeakArrayOneOneMean(kLoopnum),'-.r'); hold on
    %         plot(VTimeAWG/1e-6,WaveformIdealPeak*VIntAmpPeakArrayOneOneMean(kLoopnum),'--k','linewidth',2); hold on
    %     end
    %     plot(VTimeAWG/1e-6,VIntAmpPeakArrayOn,'b');
    %     xlabel('Time/us');
    %     ylabel('Voltage/mV');
    %     title('Envelope of output waveforms');
        


    %     subplot(2, 3, 1);
    %     if kLoopnum == NLoop
    %     % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*0.8,'-.r'); hold on
    %     plot(VTimeAWG/1e-6,WaveformIdealPeak*0.8,'--k','linewidth',2); hold on
    %     end
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm,'b');
    %     xlabel('Time/us');
    %     ylabel('Intensity/A.U.');
    %     title('Envelope of electronic waveforms');



        % Save data
        if DataSaveSaveData == 'Yes'
            save([get(handles.edit23, 'String'),'\',PresentFileName,'_AWG_pattern_Loop_',num2str(kLoopnum)],'WAVEFORMDATAARRAYNorm');
        end

        % % Process pulses that are One!!!
        % VIntAmpPeakArrayOnOne = VIntAmpPeakArray.*WavefommIdealPulseOnOne;           % Choose only area with pulse on
        % % VIntAmpPeakArrayOneOneMean(kLoopnum) = sum(VIntAmpPeakArrayOnOne)/NumPulseOnOne;   % Average value
        % VIntAmpPeakArrayOnOneNorm = VIntAmpPeakArrayOnOne/max(VIntAmpPeakArrayOnOne);     % Normalize, [0, 1]
        % VIntAmpPeakArrayOnOneNormMean = sum(VIntAmpPeakArrayOnOneNorm)/NumPulseOnOne;       % Mean value of pulse peak
        % VIntAmpPeakArrayOnOneFluc = abs(VIntAmpPeakArrayOnOneNorm-VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
        % VIntAmpPeakArrayOnOneFlucRatio = VIntAmpPeakArrayOnOneFluc/VIntAmpPeakArrayOnOneNormMean;   % Fluctuation ratio
        % VIntAmpPeakArrayOnOneAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio)/NumPulseOnOne;            % average fluctuation
        % VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne;         % RMS fluctuation
        % VIntAmpPeakArrayOnOneMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnOneFlucRatio);
        % 
        % % Process pulses that are not off
        % VIntAmpPeakArrayOn = VIntAmpPeakArray.*WavefommIdealPulseOn;           % Choose only area with pulse on
        % VIntAmpPeakArrayMeanOn(kLoopnum) = sum(VIntAmpPeakArrayOn)/NumPulseOn;   % Average value
        % VIntAmpPeakArrayOnNorm = VIntAmpPeakArrayOn/max(VIntAmpPeakArrayOn);     % Normalize, [0, 1]
        % VIntAmpPeakArrayOnNormMean = sum(VIntAmpPeakArrayOnNorm)/NumPulseOn;       % Mean value of pulse peak
        % VIntAmpPeakArrayOnFluc = abs(VIntAmpPeakArrayOnNorm-VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak); % Deviation from ideal peak
        % 
        % VIntAmpPeakArrayOnFlucRatio = VIntAmpPeakArrayOnFluc/VIntAmpPeakArrayOnNormMean;   % Fluctuation ratio
        % VIntAmpPeakArrayOnAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio)/NumPulseOn;            % average fluctuation
        % VIntAmpPeakArrayOnRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio.^2)/NumPulseOn;         % RMS fluctuation
        % VIntAmpPeakArrayOnMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnFlucRatio);

        % % use measured waveform to generate a new waveform for AWG input: WAVEFORMDATAARRAY.
        WaveformComp = VIntAmpPeakArrayOnOneNorm - VIntAmpPeakArrayOnOneNormMean.*WaveformIdealPeak;
        WaveformComp = WaveformComp.*WavefommIdealPulseOnOne;
        %if mod(kLoopnum,2) == 1 
        %    figure; plot(VTimeAWG, WaveformComp)
        %end
        %disp(VIntAmpPeakArrayOnOneNormMean)
        % use measured waveform to generate a new waveform for AWG input: WAVEFORMDATAARRAY.
%         WaveformComp = VIntAmpPeakArrayOnNorm - VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak;
%         WaveformComp = WaveformComp.*WavefommIdealPulseOn;

        % figure
        % plot(VTimeAWG,VIntAmpPeakArrayOnNorm,'b');hold on
        % plot(VTimeAWG,VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak,'k');hold on
        % plot(VTimeAWG,WaveformComp,'r');hold on
        % 
        % figure (11)
        
        
        
%         WaveformCompFactor = 4;
%         if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=0.5/100
%             WaveformCompFactor = 8;
%         end
%         if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=0.2/100
%             WaveformCompFactor = 16;
%         end
%         if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=0.02/100
%             WaveformCompFactor = 32;
%         end
%         
% 
% 
%         WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm - WaveformComp/WaveformCompFactor;
%         
%         % In case: minimum of the waveform is less than 0.
%         if min(WAVEFORMDATAARRAYNorm)<0
%             WAVEFORMDATAARRAYNorm = (WAVEFORMDATAARRAYNorm-min(WAVEFORMDATAARRAYNorm)) / max((WAVEFORMDATAARRAYNorm-min(WAVEFORMDATAARRAYNorm)));
%         else
%             WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm / max(WAVEFORMDATAARRAYNorm);
%         end
        
        
                WaveformCompFactorArray = [
        1           6
        24/100      8
        4/100       12
        2/100       16
        1/100       16
%         1           6
%         25/100       8
%         16/100       16
%         10/100       32
%         1           6
%         6/100       8
%         2/100       16
%         1/100       32
%         1           24
%         25/100       28
%         10/100     16
%         1/100    32

        ];
        % figure (11)
        WaveformCompFactor = WaveformCompFactorArray(1,2);
        if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=WaveformCompFactorArray(2,1)
            WaveformCompFactor = WaveformCompFactorArray(2,2);
        end
        if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=WaveformCompFactorArray(3,1)
            WaveformCompFactor = WaveformCompFactorArray(3,2);
        end
        if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=WaveformCompFactorArray(4,1)
            WaveformCompFactor = WaveformCompFactorArray(4,2);
        end
        if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)<=WaveformCompFactorArray(5,1)
            WaveformCompFactor = WaveformCompFactorArray(5,2);
        end
        
% CodeTest 2
%         WaveformComp = WaveformComp./max(WaveformComp);

%         % previous code from Bo 20190622
%         WAVEFORMDATAARRAYNorm  = WAVEFORMDATAARRAYNorm - WaveformComp/WaveformCompFactor;
%         WAVEFORMDATAARRAYNorm(WAVEFORMDATAARRAYNorm<0)=0;
%         WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm / max(WAVEFORMDATAARRAYNorm);
        
      % Test 20220524, AWG range 1-6 V. Does not work. Keep AWG_floor equal
      % to 0.0.
        AWG_floor = 0.0;
        AWG_floor_array = AWG_floor.*WavefommIdealPulseOnOne;
        WAVEFORMDATAARRAYNorm  = WAVEFORMDATAARRAYNorm.* (1 - WaveformComp/WaveformCompFactor);  
        WAVEFORMDATAARRAYNorm(WAVEFORMDATAARRAYNorm<0)=0;
        WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm / max(WAVEFORMDATAARRAYNorm);
        WAVEFORMDATAARRAYNorm = AWG_floor_array + (1-AWG_floor).*WAVEFORMDATAARRAYNorm;
        fprintf(['WAVEFORMDATAARRAYNorm max = ',num2str(round(max(WAVEFORMDATAARRAYNorm*1000))/1000),' \n  ']);
        fprintf(['WAVEFORMDATAARRAYNorm min = ',num2str(round(min(WAVEFORMDATAARRAYNorm*1000))/1000),' \n  ']);
        
        
        % WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm*2-1;                  % Generate the signal for AWG
        % WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm;                  % Generate the signal for AWG
        % WAVEFORMDATAARRAYNorm(WAVEFORMDATAARRAYNorm<0)=0;
        % WAVEFORMDATAARRAYNorm = abs((WAVEFORMDATAARRAY-PulseOff)/max(WAVEFORMDATAARRAY-PulseOff)); % Normalized
        
        % xlim([100,106]);
        % Share final DataArray
        handles.LastWAVEFORMDATAARRAY     = WAVEFORMDATAARRAY;
        FinalPatternInSampleFirst = PatternInSampleFirst;
        FinalPatternInSampleLast  = PatternInSampleLast;
        % First part, line 2 to 511
        for k = 1:LineOfThePatternFirst
            StartL = 1+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
            EndL   = SamplesOfOneLineFirst+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
            FinalPatternInSampleFirst(k,:)=WAVEFORMDATAARRAY(StartL:EndL);
        end
        % Second part, line 512 to 1
        StartL = 1+SamplesOfOnePatternFirst;
        EndL   = SamplesOfOneLineLast+SamplesOfOnePatternFirst;
        FinalPatternInSampleLast(:)=WAVEFORMDATAARRAY(StartL:EndL);
        handles.FinalPatternInSampleFirst = FinalPatternInSampleFirst;
        handles.FinalPatternInSampleLast  = FinalPatternInSampleLast;

        WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm*(PulseOn-PulseOff)+PulseOff;
%         PulseOn = 1;                % Waveform value when pulse is on
%         PulseOff = 0;               % Waveform value when pulse is off
%         VLoops2 = linspace(0,NLoop,NLoop*1000);
        fprintf(['Processing: ',num2str(round(toc(PostProcessTakesTime)*1000)),' ms.   ']);
        
               
        % plotting envelops
        FigureTakesTime = tic;
        FigureIfShow = get(handles.text144, 'String');
        if FigureIfShow == 'Full  '
            if videoRecord == 0
                
                fh = figure (211);
                fh.Units = 'normalized'; fh.Position = [0 0 1 1];          % maximize figure 11 when run 
                subplot(2, 2, 1);
                if kLoopnum == 1
                    % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
                    plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYImageArea(FigureShowStart:FigureShowEnd)*0.8,'-.r'); 
                    hold on
                    plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WaveformIdealPeak(FigureShowStart:FigureShowEnd)*0.8,'--k','linewidth',2); 
                    hold on
                end
    %             if kLoopnum == 1 || mod(kLoopnum,20) == 0
    %                 Colork = 1-round(kLoopnum/20)/round(NLoop/20);
    %                 plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYNorm(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1]);hold on
    %             end 
                Colork = 1-(kLoopnum)/NLoop;
                plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYNorm(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1]);hold on
    %            plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYNorm(FigureShowStart:FigureShowEnd),'g');hold on
                xlabel('Time/us');
                ylabel('Intensity/A.U.');
                title('Envelope of electronic waveforms');

                subplot(2, 2, 3);
                if kLoopnum == 1
                    % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
                    plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYImageArea(FigureShowStart:FigureShowEnd)*VIntAmpPeakArrayOneOneMean(kLoopnum),'-.r'); hold on
                    plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WaveformIdealPeak(FigureShowStart:FigureShowEnd)*VIntAmpPeakArrayOneOneMean(kLoopnum),'--k','linewidth',2); hold on
                end
    %             if kLoopnum == 1 || mod(kLoopnum,20) == 0
    %                 Colork = 1-round(kLoopnum/20)/round(NLoop/20);
    %                 plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,VIntAmpPeakArray(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1]);hold on
    %             end 
                Colork = 1-(kLoopnum)/NLoop;
                plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,VIntAmpPeakArray(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1]);hold on
                plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,- WaveformComp(FigureShowStart:FigureShowEnd)*100,'color',[Colork 1 Colork]);hold on
                xlabel('Time/us');
                title('Envelope of output waveforms');

                subplot(2, 4, 3);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneMaxFluc(1:kLoopnum)*100,'-sk');hold on;
                if kLoopnum==1
                    text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneMaxFluc(1)*100),'%']);
                    xlabel('Loops');
                    ylabel('Fluctuation of signal/%');
                end
                title(['Max, ',num2str(VIntAmpPeakArrayOnOneMaxFluc(kLoopnum)*100),'%, Compensation Factor = ',num2str(WaveformCompFactor)]);

                subplot(2, 4, 4);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneRMSFluc(1:kLoopnum)*100,'-sb');hold on;
                if kLoopnum==1
                    text(0.5,max(VIntAmpPeakArrayOnOneRMSFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneRMSFluc(1)*100),'%']);
                    xlabel('Loops');
                    ylabel('Fluctuation of signal/%');
                end
                title(['RMS, ',num2str(VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)*100),'%, Compensation Factor = ',num2str(WaveformCompFactor)]);

                subplot(2, 2, 4);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOneOneMean(1:kLoopnum),'-ok');hold on;
                xlabel('Loops');
                ylabel('Mean value of output pulse peaks/a.u.');
                title('Loops vs output pulse peaks mean value');
                fprintf(['Figure: ',num2str(round(toc(FigureTakesTime)*1000)),' ms.   ']);
    %             fprintf([num2str(kLoopnum),'th rounds: ',num2str(round(toc(AWGTakesTime)*1000)),' ms.\n']);
                pause(str2double(get(handles.edit31, 'String')));
                
                % keyboard interrupt feedback loop
                isKeyPressed = ~isempty(get(fh,'CurrentCharacter'));
                if isKeyPressed
                    break
                end
                
                
            else  % if videoRecord ~= 0, i.e. if need to record a video
                fh = figure (211);
                fh.Units = 'normalized'; fh.Position = [0 0 1 1];          % maximize figure 11 when run 
                subplot(2, 2, 1);
                if kLoopnum == 1
                    % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
%                     plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYImageArea(FigureShowStart:FigureShowEnd)*0.8,'-.r'); 
%                     hold on
%                     plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WaveformIdealPeak(FigureShowStart:FigureShowEnd)*0.8,'--k','linewidth',2); 
%                     hold on
                end

                Colork = 0;
                plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,VIntModPeakArray(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1],'linewidth',2);
%                 plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYNorm(FigureShowStart:FigureShowEnd),'color',[Colork Colork 1],'linewidth',2);
%                 hold on
                xlabel('Time/us');
                ylabel('Intensity/A.U.');
                xlim([20 45])
                ylim([0 120])
                title('Envelope of optical waveforms after I.M.');
%                 title('Envelope of electronic waveforms');

                subplot(2, 2, 3);
%                 if kLoopnum == 1
%                     % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
%                     plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WAVEFORMDATAARRAYImageArea(FigureShowStart:FigureShowEnd)*VIntAmpPeakArrayOneOneMean(kLoopnum),'-.r'); hold on
%                     plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,WaveformIdealPeak(FigureShowStart:FigureShowEnd)*VIntAmpPeakArrayOneOneMean(kLoopnum),'--k','linewidth',2); hold on
%                 end

%                 Colork = 1-(kLoopnum)/NLoop;
                plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,VIntAmpPeakArray(FigureShowStart:FigureShowEnd),'color',[204/255 0 0],'linewidth',2);
%                 hold on
%                 plot(VTimeAWG(FigureShowStart:FigureShowEnd)/1e-6,- WaveformComp(FigureShowStart:FigureShowEnd)*100,'color',[Colork 1 Colork]);hold on
                xlabel('Time/us');
                ylim([0 400])
                xlim([20 45])
                title('Envelope of output waveforms');

                subplot(2, 4, 3);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneMaxFluc(1:kLoopnum)*100,'-sk');hold on;
                if kLoopnum==1
                    text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneMaxFluc(1)*100),'%']);
                    xlabel('Loops');
                    ylabel('Fluctuation of signal/%');
                end
                title(['Max, ',num2str(VIntAmpPeakArrayOnOneMaxFluc(kLoopnum)*100),'%, Compensation Factor = ',num2str(WaveformCompFactor)]);

                subplot(2, 4, 4);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneRMSFluc(1:kLoopnum)*100,'-sb');hold on;
                if kLoopnum==1
                    text(0.5,max(VIntAmpPeakArrayOnOneRMSFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneRMSFluc(1)*100),'%']);
                    xlabel('Loops');
                    ylabel('Fluctuation of signal/%');
                end
                title(['RMS, ',num2str(VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)*100),'%, Compensation Factor = ',num2str(WaveformCompFactor)]);

                subplot(2, 2, 4);
                plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOneOneMean(1:kLoopnum),'-ok');hold on;
                xlabel('Loops');
                ylabel('Mean value of output pulse peaks/a.u.');
                title('Loops vs output pulse peaks mean value');
                fprintf(['Figure: ',num2str(round(toc(FigureTakesTime)*1000)),' ms.   ']);
    %             fprintf([num2str(kLoopnum),'th rounds: ',num2str(round(toc(AWGTakesTime)*1000)),' ms.\n']);
                pause(str2double(get(handles.edit31, 'String')));
                
                % keyboard interrupt feedback loop
                isKeyPressed = ~isempty(get(fh,'CurrentCharacter'));
                if isKeyPressed
                    break
                end
                
                % capture each frame to create a movie after the loop, in the
                % end of this callback function
                F(kLoopnum) = getframe(gcf) ;
                drawnow
            end % if videoRecord == 0
            
            
        elseif FigureIfShow == 'Simple'
            FigureTakesTime = tic;
            fh = figure (211)
            subplot(1, 2, 1);
            plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneMaxFluc(1:kLoopnum)*100,'-sk');hold on;
            if kLoopnum==1
                text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneMaxFluc(1)*100),'%']);
                xlabel('Loops');
                ylabel('Fluctuation of signal/%');
            end
            title(['Max, ',num2str(VIntAmpPeakArrayOnOneMaxFluc(kLoopnum)*100),'%']);
            subplot(1, 2, 2);
            plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOneOneMean(1:kLoopnum),'-ok');hold on;
            xlabel('Loops');
            ylabel('Mean value of output pulse peaks/a.u.');
            title('Loops vs output pulse peaks mean value');
            fprintf(['Figure: ',num2str(round(toc(FigureTakesTime)*1000)),' ms.   ']);
%             fprintf([num2str(kLoopnum),'th rounds: ',num2str(round(toc(AWGTakesTime)*1000)),' ms.\n']);
            pause(str2double(get(handles.edit31, 'String')));
            
            % keyboard interrupt feedback loop
            isKeyPressed = ~isempty(get(fh,'CurrentCharacter'));
            if isKeyPressed
                break
            end
                
        end           
        fprintf([num2str(kLoopnum),'th rounds: ',num2str(round(toc(AWGTakesTime)*1000)),' ms.\n']);
        
        % if RMS reduced to a desired value, break out of the feedback loop
        % right away 20210301
        if VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) < 0.01 || VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) > 1
            break
        end
        
        
        
%     else          % used to match the end in line 2901
    %     VItensityAmp = VItensityAmp/AverageTime;
    %     VItensityMod = VItensityMod/AverageTime;
    %     % Save data
    %     if DataSaveSaveData == 'Yes'
    %         save([get(handles.edit23, 'String'),'\',PresentFileName,'_Mod_signal_Loop_',num2str(kLoopnum)],'VItensityMod');
    %         save([get(handles.edit23, 'String'),'\',PresentFileName,'_Amp_signal_Loop_',num2str(kLoopnum)],'VItensityAmp');
    %     end
    %     if FeedbackLoopRawData == 'Yes'
    %         figure (12)
    %         subplot(2,1,1);
    %         plot(VTimeOSC/1e-6,VItensityMod);
    %         % hold on;
    %         xlabel('Time/us');
    %         ylabel('Voltage/mV');
    %         xlim([0,max(VTimeOSC/1e-6)-0.0001]);
    %         subplot(2,1,2);
    %         plot(VTimeOSC/1e-6,VItensityAmp);
    %         % hold on;
    %         xlabel('Time/us');
    %         ylabel('Voltage/mV');
    %         xlim([0,max(VTimeOSC/1e-6)-0.0001]);
    %         % ylim([-0.05,1]);
    %         title('Channel A, amplified signal');
    %     end
    %     
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%%%%%%%%%%%%%%%%%%% processing waveform %%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % Here assuming AWG sampling rate 4 MHz, laser rep rate 4 MHz, OSC sampling rate 500 MHz.
    %     % AWG sampling rate is not accurately 4 MHz, so it is scaled here.
    % 
    %     % *(should be/measure)
    %     % SamRatioOscToAWGScaled =(SamRatioOscToAWG)*(15.586/15.594)*(3.12313/3.12308)*(1.2493735/1.2493756)*(3.1234383/3.1234385)*(3.1234381/3.1234383);
    %     SamRatioOscToAWGScaled = SamRatioOscToAWG;
    %     % SamRatioOscToAWGScaled =(SamRatioOscToAWGScaled)*0.999484264671746*1.007640000000000e+02/1.007580000000000e+02;
    %     % SamRatioOscToAWGScaled = SamRatioOscToAWGScaled*0.999979344326761*1.000000600804944;
    %     % SamRatioOscToAWGScaled = SamRatioOscToAWGScaled/3.1233765*3.1233756*0.999998769858928*0.999999776823323;
    %     % SamRatioOscToAWGScaled = SamRatioOscToAWGScaled/3.1233756*3.123375*0.999999838665790*1.000005575296082;
    %     % SamRatioOscToAWGScaled = SamRatioOscToAWGScaled*0.999999169444091*1.000002746970673*0.999999619555013;
    %     % SamRatioOscToAWGScaled = SamRatioOscToAWGScaled*0.999999808838426*0.999999871933274*1.000000832433816;
    %     % 4 MHz
    %     % SamRatioOscToAWGScaled
    %     % =(SamRatioOscToAWG)/200.1*200/6000.75*6000.69/(1.249338/1.2493308); % for
    %     % 4 MHz
    % 
    % 
    %     % Set the size of the waveform (equal to AWG waveform after scaling)
    %     WaveformOSCSize = round(WAVEFORMSIZE*SamRatioOscToAWGScaled);                
    %     SamRatioOscToAWGScaledInt = floor(SamRatioOscToAWGScaled/2)*2+2;
    % 
    %     % Find when does the waveform start
    %     ChannelAStartPoint = str2double(get(handles.edit29, 'String'))*1e-6;
    %     WaveChannelAStartPoint = find(VTimeOSC>=ChannelAStartPoint & VTimeOSC<(ChannelAStartPoint)+0.01e-6); % Find when does the waveform start
    % 
    % 
    % 
    %     if WaveChannelAStartPoint(1)+WaveformOSCSize-1>OSCNlengthPerChannel
    %         errordlg('OSC record length is less than the AWG waveform period','GUI Error');
    %         return;
    %     end
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    %     % Channel B process
    %     if kLoopnum==1
    %         ChannelBStartPoint = str2double(get(handles.edit30, 'String'))*1e-6;
    %         WaveChannelBStartPoint = find(VTimeOSC>ChannelBStartPoint & VTimeOSC<(ChannelBStartPoint)+0.01e-6); % Find when does the waveform start
    % 
    %         % Process the use measured waveform to measure the peaks, average amplitude of each pulse.
    %         VIntModMatrix = zeros(SamRatioOscToAWGScaledInt,WAVEFORMSIZE);
    % 
    %         % Define first value
    %         kArraytoMatrix =1;
    %         PulseStart = round(SamRatioOscToAWGScaled*(kArraytoMatrix-1))+WaveChannelBStartPoint(1);
    %         PulseEnd   = PulseStart + SamRatioOscToAWGScaledInt-1;
    %         VIntModMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix)=VItensityMod(PulseStart:PulseEnd);
    %         PulseCenter = round((PulseStart+PulseEnd)/2);
    %         % Define rest values
    %         LastOneNumber = 0;
    %         LastPulseCenter=0;
    %         for kArraytoMatrix = 2:WAVEFORMSIZE % WAVEFORMSIZE is the AWG length
    %             if WaveformIdealPeak(kArraytoMatrix-1)==1 % When last pulse is on, i.e., value is 1, then update pulse center at first
    %                 PulseCenterInt = round(PulseCenter); % Pulse center may not be int value, so round it 
    %                 LastPulseStart = PulseCenterInt-SamRatioOscToAWGScaledInt/2; % When does last pulse start, pulse center is the center of last pulse
    %                 LastPulseEnd   = PulseCenterInt+SamRatioOscToAWGScaledInt/2 -1; % When does last pulse end, pulse center is the center of last pulse
    %                 [~,PulseCenter]=max(  VItensityMod(LastPulseStart:LastPulseEnd)  ); % Find the pulse center of the last pulse, update the value with pulse peak
    %                 PulseCenter = PulseCenter + LastPulseStart-1;
    %                 PulseStart = (PulseCenter+SamRatioOscToAWGScaledInt/2);   % When does present pulse start, pulse center is the center of last pulse
    %                 PulseEnd   = (PulseCenter+SamRatioOscToAWGScaledInt*3/2) -1; % When does present pulse end, pulse center is the center of last pulse
    %                 VIntModMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix) = VItensityMod(PulseStart:PulseEnd); % Define the new pulse
    %     %             if LastOneNumber>0
    %     %                 ThisOneNumber = kArraytoMatrix-1;
    %     % %                 SamRatioOscToAWGScaled
    %     % %                 SamRatioOscToAWGScaled2 = abs(PulseCenter-0)/(ThisOneNumber-0)
    %     %             end
    %                 LastPulseCenter = PulseCenter;
    %                 LastOneNumber = kArraytoMatrix-1;
    % 
    %                 PulseCenter= PulseCenter+SamRatioOscToAWGScaled;       % update the pulse center, use the present pulse center as last pulse center for next pulse;
    %             else % When pulse is off, i.e., value is 0, then no need to update pulse center
    %                 PulseCenterInt = round(PulseCenter);    % Pulse center may not be int value, so round it 
    %                 PulseStart = (PulseCenterInt+SamRatioOscToAWGScaledInt/2);   % When does present pulse start, pulse center is the center of last pulse
    %                 PulseEnd   = (PulseCenterInt+SamRatioOscToAWGScaledInt*3/2) -1; % When does present pulse end, pulse center is the center of last pulse
    %                 VIntModMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix) = VItensityMod(PulseStart:PulseEnd); % Define the new pulse
    %                 PulseCenter= PulseCenter+SamRatioOscToAWGScaled;       % update the pulse center, use the present pulse center as last pulse center for next pulse;
    %             end
    %         end
    %         VIntModPeakArray = max(VIntModMatrix);                          % Find the peak of each pulse
    %         % Process Only One fluc
    %         VIntModPeakArrayOnOne = VIntModPeakArray.*WavefommIdealPulseOnOne;           % Choose only area with pulse on
    %         % VIntModPeakArrayOneOneMean(kLoopnum) = sum(VIntModPeakArrayOnOne)/NumPulseOnOne;   % Average value
    %         VIntModPeakArrayOnOneNorm = VIntModPeakArrayOnOne/max(VIntModPeakArrayOnOne);     % Normalize, [0, 1]
    %         VIntModPeakArrayOnOneNormMean = sum(VIntModPeakArrayOnOneNorm)/NumPulseOnOne;       % Mean value of pulse peak
    %         VIntModPeakArrayOnOneFluc = abs(VIntModPeakArrayOnOneNorm-VIntModPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
    %         VIntModPeakArrayOnOneFlucRatio = VIntModPeakArrayOnOneFluc/VIntModPeakArrayOnOneNormMean;   % Fluctuation ratio
    %         % VIntModPeakArrayOnOneAveFluc(kLoopnum) = sum(VIntModPeakArrayOnOneFlucRatio)/NumPulseOnOne;            % average fluctuation
    %         VIntModPeakArrayOnOneRMSFluc = sum(VIntModPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne;         % RMS fluctuation
    %         VIntModPeakArrayOnOneMaxFluc = max(VIntModPeakArrayOnOneFlucRatio);
    %     end
    % 
    % 
    % 
    %     % SamRatioOscToAWGScaled0 =(1000/4)/200.1*200/6000.75*6000.69;
    % 
    % 
    % 
    % 
    %     % Process the use measured waveform to measure the peaks, average amplitude of each pulse.
    %     VIntAmpMatrix = zeros(SamRatioOscToAWGScaledInt,WAVEFORMSIZE);
    %     % maxofKpulse = linspace(0,0,WAVEFORMSIZE);
    %     % Start of each pulse
    % 
    %     % WaveformIdealPeak, is AWG ideal peak
    % 
    %     % fprintf(['Test 1 takes',num2str(toc(FeedbackLoopTimeBeforeLoop)),' seconds.\n']);
    % 
    %     % WaveChannelAStartPoint = find(VTimeOSC>=ChannelAStartPoint & VTimeOSC<(ChannelAStartPoint)+0.01e-6); % Find when does the waveform start
    % 
    % 
    %     % Define first value
    %     kArraytoMatrix =1;
    %     PulseStart = round(SamRatioOscToAWGScaled*(kArraytoMatrix-1))+WaveChannelAStartPoint(1);
    %     PulseEnd   = PulseStart + SamRatioOscToAWGScaledInt-1;
    %     VIntAmpMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix)=VItensityAmp(PulseStart:PulseEnd);
    %     PulseCenter = round((PulseStart+PulseEnd)/2);
    %     % Define rest values
    %         if PulsePeakTest~=0
    %             if kArraytoMatrix ==PulsePeakTest
    %                 figure (13)
    %                 plot(VTimeOSC/1e-6,VItensityAmp);hold on;
    %                 NCheck = PulseStart:PulseEnd;
    %                 plot(VTimeOSC(NCheck)/1e-6,VItensityAmp(NCheck),'o');hold on
    %                 minVTime = min(VTimeOSC(NCheck)/1e-6)-0.1;
    %                 maxVTime = max(VTimeOSC(NCheck)/1e-6)+0.1;
    %                 xlim([minVTime,maxVTime])
    %             end
    %         end
    % 
    %     for kArraytoMatrix = 2:WAVEFORMSIZE % WAVEFORMSIZE is the AWG length
    %         if WaveformIdealPeak(kArraytoMatrix-1)==0.99 % When last pulse is on, i.e., value is 1, then update pulse center at first
    %             PulseCenterInt = round(PulseCenter); % Pulse center may not be int value, so round it 
    %             LastPulseStart = PulseCenterInt-SamRatioOscToAWGScaledInt/2; % When does last pulse start, pulse center is the center of last pulse
    %             LastPulseEnd   = PulseCenterInt+SamRatioOscToAWGScaledInt/2 -1; % When does last pulse end, pulse center is the center of last pulse
    %             [~,PulseCenter]=max(  VItensityAmp(LastPulseStart:LastPulseEnd)  ); % Find the pulse center of the last pulse, update the value with pulse peak
    %             PulseCenter = PulseCenter + LastPulseStart-1;
    %             PulseStart = (PulseCenter+SamRatioOscToAWGScaledInt/2);   % When does present pulse start, pulse center is the center of last pulse
    %             PulseEnd   = (PulseCenter+SamRatioOscToAWGScaledInt*3/2) -1; % When does present pulse end, pulse center is the center of last pulse
    %             VIntAmpMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix) = VItensityAmp(PulseStart:PulseEnd); % Define the new pulse
    %             PulseCenter= PulseCenter+SamRatioOscToAWGScaled;       % update the pulse center, use the present pulse center as last pulse center for next pulse;
    %         else % When pulse is off, i.e., value is 0, then no need to update pulse center
    %             PulseCenterInt = round(PulseCenter);    % Pulse center may not be int value, so round it 
    %             PulseStart = (PulseCenterInt+SamRatioOscToAWGScaledInt/2);   % When does present pulse start, pulse center is the center of last pulse
    %             PulseEnd   = (PulseCenterInt+SamRatioOscToAWGScaledInt*3/2) -1; % When does present pulse end, pulse center is the center of last pulse
    %             VIntAmpMatrix(1:SamRatioOscToAWGScaledInt,kArraytoMatrix) = VItensityAmp(PulseStart:PulseEnd); % Define the new pulse
    %             PulseCenter= PulseCenter+SamRatioOscToAWGScaled;       % update the pulse center, use the present pulse center as last pulse center for next pulse;
    %         end
    %     % Use this to test measuring of pulse peak
    %     % max = 1000000;
    %         if PulsePeakTest~=0
    %             if kArraytoMatrix == PulsePeakTest
    %                 figure (13)
    %                 plot(VTimeOSC/1e-6,VItensityAmp);hold on;
    %                 NCheck = PulseStart:PulseEnd;
    % 
    %                 plot(VTimeOSC(NCheck)/1e-6,VItensityAmp(NCheck),'o');hold on
    %                 minVTime = min(VTimeOSC(NCheck)/1e-6)-0.1;
    %                 maxVTime = max(VTimeOSC(NCheck)/1e-6)+0.1;
    %                 xlim([minVTime,maxVTime])
    %             end
    %         end
    % 
    % 
    %     end
    %     % fprintf(['Test 1 takes',num2str(toc(FeedbackLoopTimeBeforeLoop)),' seconds.\n']);
    %     % plot(VTimeAWG/1e-6,maxofKpulse);hold on
    % 
    %     % xlabel('Time/us');
    % 
    % 
    %     % whos VIntAmpMatrix
    % 
    %     % whos VIntAmpMatrix
    % 
    %     VIntAmpPeakArray = max(VIntAmpMatrix);                          % Find the peak of each pulse
    %     % whos VIntAmpPeakArray
    %     % whos VIntAmpPeakArray
    % 
    %     % Process pulses that are One!!!
    %     VIntAmpPeakArrayOnOne = VIntAmpPeakArray.*WavefommIdealPulseOnOne;           % Choose only area with pulse on
    %     VIntAmpPeakArrayOneOneMean(kLoopnum) = sum(VIntAmpPeakArrayOnOne)/NumPulseOnOne;   % Average value
    %     VIntAmpPeakArrayOnOneNorm = VIntAmpPeakArrayOnOne/max(VIntAmpPeakArrayOnOne);     % Normalize, [0, 1]
    %     VIntAmpPeakArrayOnOneNormMean = sum(VIntAmpPeakArrayOnOneNorm)/NumPulseOnOne;       % Mean value of pulse peak
    %     VIntAmpPeakArrayOnOneFluc = abs(VIntAmpPeakArrayOnOneNorm-VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
    %     VIntAmpPeakArrayOnOneFlucRatio = VIntAmpPeakArrayOnOneFluc/VIntAmpPeakArrayOnOneNormMean;   % Fluctuation ratio
    %     VIntAmpPeakArrayOnOneAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio)/NumPulseOnOne;            % average fluctuation
    %     VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne;         % RMS fluctuation
    %     VIntAmpPeakArrayOnOneMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnOneFlucRatio);
    % 
    %     % Process pulses that are not off
    %     VIntAmpPeakArrayOn = VIntAmpPeakArray.*WavefommIdealPulseOn;           % Choose only area with pulse on
    %     VIntAmpPeakArrayMeanOn(kLoopnum) = sum(VIntAmpPeakArrayOn)/NumPulseOn;   % Average value
    %     VIntAmpPeakArrayOnNorm = VIntAmpPeakArrayOn/max(VIntAmpPeakArrayOn);     % Normalize, [0, 1]
    %     VIntAmpPeakArrayOnNormMean = sum(VIntAmpPeakArrayOnNorm)/NumPulseOn;       % Mean value of pulse peak
    %     VIntAmpPeakArrayOnFluc = abs(VIntAmpPeakArrayOnNorm-VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak); % Deviation from ideal peak
    % 
    %     VIntAmpPeakArrayOnFlucRatio = VIntAmpPeakArrayOnFluc/VIntAmpPeakArrayOnNormMean;   % Fluctuation ratio
    %     VIntAmpPeakArrayOnAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio)/NumPulseOn;            % average fluctuation
    %     VIntAmpPeakArrayOnRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio.^2)/NumPulseOn;         % RMS fluctuation
    %     VIntAmpPeakArrayOnMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnFlucRatio);
    % 
    %     % figure
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnNorm,'b');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak,'k');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnFluc,'r');hold on
    %     % plot(VTimeAWG,WaveformIdealPeak,'--g');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne,'k');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnOneFluc,'r');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOn,'--k');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnOneNorm,'--g');hold on
    % 
    %     % figure
    %     % plot(VTimeOSC/1e-6-ChannelAStartPoint/1e-6,VItensityAmp);hold on;
    %     % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm*20,'k');hold on
    % 
    %     % hold on;
    %     % xlabel('Time/us');
    %     % ylabel('Voltage/mV');
    %     % 
    %     % fprintf(['Test 2 takes',num2str(toc(FeedbackLoopTimeBeforeLoop)),' seconds.\n']);
    %     % % 
    %     % figure (13)
    %     % plot(VTimeAWG/1e-6,VTimeAWG.*0+VIntAmpPeakArrayOnNormMean,'--r');hold on
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnNorm,'-.k');hold on
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnFluc,'-b');hold on
    %     % xlim([0,40]);
    %     figure (11)
    %     % set(gcf, 'units','normalized','outerposition',[0 0.1 0.5 0.9]);
    %     % Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
    %     % caption = sprintf('Original "coins" image showing\n6 nickels (the larger coins) and 4 dimes (the smaller coins).');
    %     % title(caption, 'FontSize', 14);
    %     % axis image; % Make sure image is not artificially stretched because of screen's
    %     % 
    %     % % figure (10)
    %     subplot(2, 3, 4);
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnNorm);hold on
    %     if kLoopnum == NLoop
    %         % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %         plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*VIntAmpPeakArrayOneOneMean(kLoopnum),'-.r'); hold on
    %         plot(VTimeAWG/1e-6,WaveformIdealPeak*VIntAmpPeakArrayOneOneMean(kLoopnum),'--k','linewidth',2); hold on
    %     end
    %     plot(VTimeAWG/1e-6,VIntAmpPeakArrayOn,'b');
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOn,'--k');hold on
    %     % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm*40,'-.r');hold on
    %     xlabel('Time/us');
    %     % ylabel('Intensity/A.U.');
    %     title('Envelope of output waveforms');
    %     % legend('1st loop','2nd','3rd','4th','...');
    %     % xlim([960,1010]);
    %     % xlim([10,50]);
    %     % % Save data
    % 
    % 
    % 
    %     % xlim([7570,7630]);
    %     % figure (10)
    % 
    %     % subplot(4, 3, 4);
    %     % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm);
    %     % xlabel('Time/us');
    %     % ylabel('AWG waveform');
    %     % ylim([-0.05,1.05]);
    %     % % title('Envelope of electronic waveforms');
    %     % 
    %     % subplot(4, 3, 1);
    %     % plot(VTimeAWG/1e-6,VIntModPeakArray);
    %     % 
    %     % % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOn,'--k');hold on
    %     % % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm*40,'-.r');hold on
    %     % xlabel('Time/us');
    %     % ylabel('Initial modulated light');
    %     % ylabel('Intensity/A.U.');
    %     % title('Envelope of initial modulated waveforms');
    % 
    %     % xlim([960,1010]);
    %     % xlim([10,50]);
    %     % legend('1st loop','2nd','3rd','4th','...');
    %     
    %     
    %     if kLoopnum == NLoop
    %         % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %         plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*0.8,'-.r'); hold on
    %         plot(VTimeAWG/1e-6,WaveformIdealPeak*0.8,'--k','linewidth',2); hold on
    %     end
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm,'b');
    %     xlabel('Time/us');
    %     ylabel('Intensity/A.U.');
    %     title('Envelope of electronic waveforms');
    %     xlim([0,300]);
    %     xlim([100,106]);
    %     % xlim([7570,7630]);
    % 
    %     subplot(2, 3, 1);
    %     if kLoopnum == NLoop
    %     % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*0.8,'-.r'); hold on
    %     plot(VTimeAWG/1e-6,WaveformIdealPeak*0.8,'--k','linewidth',2); hold on
    %     end
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm,'b');
    %     xlabel('Time/us');
    %     ylabel('Intensity/A.U.');
    %     title('Envelope of electronic waveforms');
    % 
    % 
    % 
    %     % Save data
    %     if DataSaveSaveData == 'Yes'
    %         save([get(handles.edit23, 'String'),'\',PresentFileName,'_AWG_pattern_Loop_',num2str(kLoopnum)],'WAVEFORMDATAARRAYNorm');
    %     end
    % 
    % 
    %     % % Process pulses that are One!!!
    %     % VIntAmpPeakArrayOnOne = VIntAmpPeakArray.*WavefommIdealPulseOnOne;           % Choose only area with pulse on
    %     % % VIntAmpPeakArrayOneOneMean(kLoopnum) = sum(VIntAmpPeakArrayOnOne)/NumPulseOnOne;   % Average value
    %     % VIntAmpPeakArrayOnOneNorm = VIntAmpPeakArrayOnOne/max(VIntAmpPeakArrayOnOne);     % Normalize, [0, 1]
    %     % VIntAmpPeakArrayOnOneNormMean = sum(VIntAmpPeakArrayOnOneNorm)/NumPulseOnOne;       % Mean value of pulse peak
    %     % VIntAmpPeakArrayOnOneFluc = abs(VIntAmpPeakArrayOnOneNorm-VIntAmpPeakArrayOnOneNormMean.*WavefommIdealPulseOnOne); % Deviation from ideal peak
    %     % VIntAmpPeakArrayOnOneFlucRatio = VIntAmpPeakArrayOnOneFluc/VIntAmpPeakArrayOnOneNormMean;   % Fluctuation ratio
    %     % VIntAmpPeakArrayOnOneAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio)/NumPulseOnOne;            % average fluctuation
    %     % VIntAmpPeakArrayOnOneRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnOneFlucRatio.^2)/NumPulseOnOne;         % RMS fluctuation
    %     % VIntAmpPeakArrayOnOneMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnOneFlucRatio);
    %     % 
    %     % % Process pulses that are not off
    %     % VIntAmpPeakArrayOn = VIntAmpPeakArray.*WavefommIdealPulseOn;           % Choose only area with pulse on
    %     % VIntAmpPeakArrayMeanOn(kLoopnum) = sum(VIntAmpPeakArrayOn)/NumPulseOn;   % Average value
    %     % VIntAmpPeakArrayOnNorm = VIntAmpPeakArrayOn/max(VIntAmpPeakArrayOn);     % Normalize, [0, 1]
    %     % VIntAmpPeakArrayOnNormMean = sum(VIntAmpPeakArrayOnNorm)/NumPulseOn;       % Mean value of pulse peak
    %     % VIntAmpPeakArrayOnFluc = abs(VIntAmpPeakArrayOnNorm-VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak); % Deviation from ideal peak
    %     % 
    %     % VIntAmpPeakArrayOnFlucRatio = VIntAmpPeakArrayOnFluc/VIntAmpPeakArrayOnNormMean;   % Fluctuation ratio
    %     % VIntAmpPeakArrayOnAveFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio)/NumPulseOn;            % average fluctuation
    %     % VIntAmpPeakArrayOnRMSFluc(kLoopnum) = sum(VIntAmpPeakArrayOnFlucRatio.^2)/NumPulseOn;         % RMS fluctuation
    %     % VIntAmpPeakArrayOnMaxFluc(kLoopnum) = max(VIntAmpPeakArrayOnFlucRatio);
    % 
    % 
    %     % % use measured waveform to generate a new waveform for AWG input: WAVEFORMDATAARRAY.
    %     WaveformComp = VIntAmpPeakArrayOnOneNorm - VIntAmpPeakArrayOnOneNormMean.*WaveformIdealPeak;
    %     WaveformComp = WaveformComp.*WavefommIdealPulseOnOne;
    % 
    % 
    %     % use measured waveform to generate a new waveform for AWG input: WAVEFORMDATAARRAY.
    %     % WaveformComp = VIntAmpPeakArrayOnNorm - VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak;
    %     % WaveformComp = WaveformComp.*WavefommIdealPulseOn;
    % 
    %     % figure
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnNorm,'b');hold on
    %     % plot(VTimeAWG,VIntAmpPeakArrayOnNormMean.*WaveformIdealPeak,'k');hold on
    %     % plot(VTimeAWG,WaveformComp,'r');hold on
    %     % 
    %     % figure (11)
    %     WaveformCompFactor = 4;
    %     if VIntAmpPeakArrayOnRMSFluc(kLoopnum)<=0.5/100
    %         WaveformCompFactor = 8;
    %     end
    %     if VIntAmpPeakArrayOnRMSFluc(kLoopnum)<=0.05/100
    %         WaveformCompFactor = 16;
    %     end
    %     if VIntAmpPeakArrayOnRMSFluc(kLoopnum)<=0.01/100
    %         WaveformCompFactor = 32;
    %     end
    %     WaveformCompFactor = 6;
    % 
    %     WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm - WaveformComp/WaveformCompFactor;
    % 
    %     if min(WAVEFORMDATAARRAYNorm)<0
    %         WAVEFORMDATAARRAYNorm = (WAVEFORMDATAARRAYNorm-min(WAVEFORMDATAARRAYNorm)) / max((WAVEFORMDATAARRAYNorm-min(WAVEFORMDATAARRAYNorm)));
    %     else
    %         WAVEFORMDATAARRAYNorm = WAVEFORMDATAARRAYNorm / max(WAVEFORMDATAARRAYNorm);
    % 
    %     end
    %     % WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm*2-1;                  % Generate the signal for AWG
    %     % WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm;                  % Generate the signal for AWG
    %     % WAVEFORMDATAARRAYNorm(WAVEFORMDATAARRAYNorm<0)=0;
    %     % WAVEFORMDATAARRAYNorm = abs((WAVEFORMDATAARRAY-PulseOff)/max(WAVEFORMDATAARRAY-PulseOff)); % Normalized
    %     % 
    %     subplot(2, 3, 5);
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOnNorm);hold on
    % 
    %     if kLoopnum == NLoop
    %     % plot(VTimeAWG/1e-6,VIntModPeakArray,'-.k');hold on
    %     plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYImageArea*VIntAmpPeakArrayOneOneMean(kLoopnum),'-.r'); hold on
    %     plot(VTimeAWG/1e-6,WaveformIdealPeak*VIntAmpPeakArrayOneOneMean(kLoopnum),'--k','linewidth',2); hold on
    %     end
    % 
    %     plot(VTimeAWG/1e-6,VIntAmpPeakArray,'b');hold on
    %     plot(VTimeAWG/1e-6,- WaveformComp*100,'g');hold on
    %     % xlim([100,106]);
    %     % plot(VTimeAWG/1e-6,VIntAmpPeakArrayOn,'--k');hold on
    %     % plot(VTimeAWG/1e-6,WAVEFORMDATAARRAYNorm*40,'-.r');hold on
    %     xlabel('Time/us');
    %     % ylabel('Intensity/A.U.');
    %     title('Envelope of output waveforms');
    %     % xlim([0,300]);
    % 
    % 
    %     xlim([100,106]);
    %     % Share final DataArray
    %     handles.LastWAVEFORMDATAARRAY     = WAVEFORMDATAARRAY;
    % 
    %     FinalPatternInSampleFirst = PatternInSampleFirst;
    %     FinalPatternInSampleLast  = PatternInSampleLast;
    %     % First part, line 2 to 511
    %     for k = 1:LineOfThePatternFirst
    %         StartL = 1+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
    %         EndL   = SamplesOfOneLineFirst+(k-1)*(SamplesOfOneLineFirst+Scan_Comp_Sample);
    %         FinalPatternInSampleFirst(k,:)=WAVEFORMDATAARRAY(StartL:EndL);
    %     end
    %     % Second part, line 512 to 1
    %     StartL = 1+SamplesOfOnePatternFirst;
    %     EndL   = SamplesOfOneLineLast+SamplesOfOnePatternFirst;
    %     FinalPatternInSampleLast(:)=WAVEFORMDATAARRAY(StartL:EndL);
    %     handles.FinalPatternInSampleFirst = FinalPatternInSampleFirst;
    %     handles.FinalPatternInSampleLast  = FinalPatternInSampleLast;
    % 
    %     WAVEFORMDATAARRAY = WAVEFORMDATAARRAYNorm*(PulseOn-PulseOff)+PulseOff;
    % 
    % 
    %     % whos WAVEFORMDATAARRAY;
    %     % MAXW = max(WAVEFORMDATAARRAYNorm)
    %     % MINW = min(WAVEFORMDATAARRAYNorm)
    % 
    % 
    % 
    % 
    %     PulseOn = 1;                % Waveform value when pulse is on
    %     PulseOff = 0;               % Waveform value when pulse is off
    %     VLoops2 = linspace(0,NLoop,NLoop*1000);
    %     subplot(2, 6, 5);
    %     plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneMaxFluc(1:kLoopnum)*100,'-sk');hold on;
    %     if kLoopnum==1
    %         plot(VLoops2,VLoops2.*0+VIntModPeakArrayOnOneMaxFluc*100,'-.k');hold on;
    %         text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneMaxFluc(1)*100),'%']);
    %         text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*7*100,['Initial Mod = ',num2str(VIntModPeakArrayOnOneMaxFluc*100),'%']);
    %         xlabel('Loops');
    %         ylabel('Fluctuation of signal/%');
    %     %     title('Max Fluctuation');
    %     end
    %     if kLoopnum==NLoop
    %         legend('Amp Max','Mod Max');
    %     end
    %     % text(0.5,max(VIntAmpPeakArrayOnOneMaxFluc(1))/10*6*100,['Present Amp = ',num2str(VIntAmpPeakArrayOnOneMaxFluc(kLoopnum)*100),'%']);
    %     title(['Max, ',num2str(VIntAmpPeakArrayOnOneMaxFluc(kLoopnum)*100),'%']);
    % 
    %     subplot(2, 6, 6);
    %     plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOnOneRMSFluc(1:kLoopnum)*100,'-sb');hold on;
    % 
    %     if kLoopnum==1
    %         plot(VLoops2,0.*VLoops2+VIntModPeakArrayOnOneRMSFluc*100,'-.b');hold on;
    %         text(0.5,max(VIntAmpPeakArrayOnOneRMSFluc(1))/10*8*100,['Initial Amp = ',num2str(VIntAmpPeakArrayOnOneRMSFluc(1)*100),'%']);
    %         text(0.5,max(VIntAmpPeakArrayOnOneRMSFluc(1))/10*7*100,['Initial Mod = ',num2str(VIntModPeakArrayOnOneRMSFluc*100),'%']);
    %         xlabel('Loops');
    %         ylabel('Fluctuation of signal/%');
    % 
    %     end
    %     if kLoopnum==NLoop
    %         legend('Amp RMS','Mod RMS');
    %     end
    %     % text(0.5,max(VIntAmpPeakArrayOnOneRMSFluc(1))/10*6*100,['Present Amp = ',num2str(VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)*100),'%']);
    %     title(['RMS, ',num2str(VIntAmpPeakArrayOnOneRMSFluc(kLoopnum)*100),'%']);
    % 
    % 
    %     subplot(2, 3, 6);
    %     plot(VLoops(1:kLoopnum),VIntAmpPeakArrayOneOneMean(1:kLoopnum),'-ok');hold on;
    %     xlabel('Loops');
    %     ylabel('Mean value of output pulse peaks/a.u.');
    %     title('Loops vs output pulse peaks mean value');
    % 
    % 
    %     fprintf([num2str(kLoopnum),'th rounds, ']);
    %     pause(str2double(get(handles.edit31, 'String')));

%     end                   % used to match the if in line 2008 


end

% Save data
if DataSaveSaveData == 'Yes'
save([get(handles.edit23, 'String'),'\',PresentFileName,'_OSC_Time'],'VTimeOSC');
save([get(handles.edit23, 'String'),'\',PresentFileName,'_AWG_Time'],'VTimeAWG');
set(handles.text38, 'String',get(handles.text37, 'String'));
set(handles.text37, 'String',get(handles.text36, 'String'));
set(handles.text36, 'String',PresentFileName);
end

fprintf(['Feedback loop takes',num2str(toc(FeedbackLoopTime)),' seconds.\n']);
% figure (13)


% export AWG waveform that results in equalized optical output
csvwrite('patternAfterFdbk.csv', WAVEFORMDATAARRAY);

if videoRecord == 1
    % a movie capturing each loop
    % create the video writer with 2 fps
    writerObj = VideoWriter('myVideo.avi');
    writerObj.FrameRate = 2;

    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end


NewhObject=hObject;
Newhandles=handles;

% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DataSaveSavedDataRead(handles);

% 4.1 Data save
function DataSaveSavedDataRead(handles)

Folder = [get(handles.edit33, 'String')];
FileName = [get(handles.edit34, 'String')];
FinalName = [Folder,'/',FileName,'.mat'];
ReadData = importdata(FinalName);
if FileName(1:12)=='OSC_Measure_'
    VTimeOSC = ReadData.VTimeOSC;
    VItensityMod=ReadData.VItensityMod;
    VItensityAmp = ReadData.VItensityAmp;
    figure
subplot(2,1,1);
plot(VTimeOSC/1e-6,VItensityMod);
xlabel('Time/us');
ylabel('Voltage/mV');
xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
title('Channel B');

subplot(2,1,2);
plot(VTimeOSC/1e-6,VItensityAmp);
xlabel('Time/us');
ylabel('Voltage/mV');
xlim([min(VTimeOSC/1e-6),max(VTimeOSC/1e-6)-0.0001]);
% ylim([-0.05,1]);
title('Channel A');
else
    
    
    
end



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Initial AWG
[hObject,handles] = AWGInitializationFunction(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);
% AWG Initialization status update
AWGInitializationCheck(handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = ImageProcessMainFunction(handles);
% Update handles structure
guidata(hObject, handles);

function edit11_Callback(~, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% OSC Initialization status update
OSCInitializationCheck(handles);


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Initialization status update
AWGInitializationCheck(handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataSaveSaveModeNum = handles.DataSaveSaveModeNum+1;
NewNum = rem(handles.DataSaveSaveModeNum,2)+1;
handles.DataSaveSaveMode = ['Yes';'No '];
set(handles.text25, 'String', handles.DataSaveSaveMode(NewNum,:));
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.text24, 'String', 'Test mode'); 







% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('GUI is fine!');



% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles]=AWGCheckAndRemove(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);
% Initialization status update
AWGInitializationCheck(handles);




% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[hObject, handles]=OSCMeasureTest(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);
% OSC Initialization status update
OSCInitializationCheck(handles);




function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit17_Callback(hObject, ~, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, ~, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, ~, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PresentText = get(handles.text26, 'String'); 
if PresentText == 'External'
    set(handles.text26, 'String', 'Internal');
else
    set(handles.text26, 'String', 'External');
end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, ~, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, ~, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.text24, 'String', 'Work mode'); 






function edit20_Callback(hObject, ~, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, ~, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, ~, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, ~, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, ~, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, ~, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 5.1 Check if image exist
function ImageProcessImageExistence(FileName, handles)

if ~exist(FileName, 'file')
    errordlg('Loaded image file does not exist','GUI Error');
    return;
end

FileTif= FileName;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);

set(handles.text73,'string',num2str(mImage));
set(handles.text78,'string',num2str(nImage));
set(handles.text80,'string',num2str(NumberImages));
if NumberImages>3
    set(handles.text70,'string','Stack');
else
    set(handles.text70,'string','Image');
end



% % 5.2 Stack fluctuation process
function ImageProcessed = ImageProcessStackProcess(handles)

fullFileName = [get(handles.edit25,'string'),'\',get(handles.edit24,'string')];


FileTif= fullFileName;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
 
FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
for i=1:NumberImages
    TemReadImage = imread(FileTif,'Index',i);
    if size(TemReadImage,3)==3
        TemReadImage = rgb2gray(TemReadImage);
    end
    FinalImage(:,:,i) = TemReadImage;
%    FinalImage(:,:,i)=imread(FileTif,'Index',i);
end

FinalImagePro = permute(FinalImage,[3 1 2]);
FImax = max(FinalImagePro);
FImin = min(FinalImagePro);
FIFluc = FImax-FImin;
FIFluc2 = permute(FIFluc,[2 3 1]);
[pixelCount, grayLevels] = imhist(FIFluc2);
MeanPixelValue = sum(pixelCount)/max(grayLevels);
[pixelCountNumber,~]=find(pixelCount>MeanPixelValue/2);
CutValue = max(pixelCountNumber);

NormalizationFactor = str2double(get(handles.edit38,'string'));
FIFluc = FIFluc2*( max(grayLevels)/CutValue * NormalizationFactor);


% If we get here, we should have found the image file.
originalImage = FIFluc;
ImageProcessed = originalImage;

% 5.3 Single image process
function ImageProcessed = ImageProcessImageProcess(handles)

fullFileName = [get(handles.edit25,'string'),'\',get(handles.edit24,'string')];

FinalImage = imread(fullFileName);
% if rgb, then convert to gray
if size(FinalImage,3) == 3
    FinalImage = rgb2gray(FinalImage);
end

% verticalLineShift = 1;              % number of vertical line shift, due to the delay physically added in polygon loop closure
% vls = verticalLineShift;
% Image3 = zeros(size(FinalImage));
% Image3(1:end-vls,:) = FinalImage(1+vls:end,:);
% FinalImage = Image3;


[pixelCount, grayLevels] = imhist(FinalImage);
MeanPixelValue = sum(pixelCount)/max(grayLevels);
[pixelCountNumber,~]=find(pixelCount>MeanPixelValue/2);
CutValue = max(pixelCountNumber);
NormalizationFactor = str2double(get(handles.edit38,'string'));
FinalImage2 = FinalImage*( max(grayLevels)/CutValue * NormalizationFactor);


% If we get here, we should have found the image file.
originalImage = FinalImage2;
ImageProcessed = originalImage;

% 5.x main function
function Newhandles = ImageProcessMainFunction(handles)

captionFontSize = 14;
% RelationshipPixelAndSampes = importdata('Pixel_vs_samples.txt');
% RelationshipPixelAndSampes = importdata('Pixel_vs_samples20190118.txt');
RelationshipPixelAndSampes = importdata(get(handles.edit54,'string'));

% Check if the file exist
fullFileName = [get(handles.edit25,'string'),'\',get(handles.edit24,'string')];
ImageProcessImageExistence(fullFileName, handles);


% Stack or Image?
FileType = get(handles.text70,'string');
if FileType == 'Stack'
    originalImage = ImageProcessStackProcess(handles);
else
    originalImage = ImageProcessImageProcess(handles);
end

% If figure 88 exist, then close it.
if ishandle(88)
    close figure 88;
end
figure (88)
set(gcf, 'units','normalized','outerposition',[0.1 0 0.9 0.8]);
% Maximize the figure window.
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
axis image; % Make sure image is not artificially stretched because of screen's
% Display the grayscale image.
subplot(2, 3, 2);
imshow(originalImage);
caption = sprintf('Fluctuation of each pixel in the stack');
title(caption, 'FontSize', captionFontSize);

% Just for fun, let's get its histogram and display it.
[pixelCount, grayLevels] = imhist(originalImage);
figure (88)
subplot(4, 3, 1);
bar(pixelCount);
title('Histogram of original image', 'FontSize', captionFontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.
grid on;
hold on;


% Method #2: using a logical operation.
thresholdValue = str2double(get(handles.edit35,'string'));
binaryImage = originalImage > thresholdValue; % Bright objects will be chosen if you use >.
% Show the threshold as a vertical red bar on the histogram.

maxYValue = ylim;
line([thresholdValue, thresholdValue], maxYValue, 'Color', 'r');
% Place a text label on the bar chart showing the threshold.
annotationText = sprintf('Thresholded at %d gray levels', thresholdValue);
% For text(), the x and y need to be of the data class "double" so let's cast both to double.
text(double(thresholdValue + 5), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
text(double(thresholdValue - 70), double(0.94 * maxYValue(2)), 'Background', 'FontSize', 10, 'Color', [0 0 .5]);
text(double(thresholdValue + 30), double(0.94 * maxYValue(2)), 'Foreground', 'FontSize', 10, 'Color', [0 0 .5]);

% Display the binary image.
figure (88)
subplot(2, 3, 3);
imshow(binaryImage); 
title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize); 

% Remove small features, Display the processed binary image.
FilterPixelForBinaryImage = str2double(get(handles.edit37,'string'));
binaryImage = bwareaopen(binaryImage,FilterPixelForBinaryImage);
figure (88)
subplot(2, 3, 6);
imshow(binaryImage); 
title('Binary Image, removed small features', 'FontSize', captionFontSize); 

% Image grows, in case of mouse motion.
ImageGrowPixelsNumber = str2double(get(handles.edit36,'string'));
for k=1:ImageGrowPixelsNumber
BW=edge(binaryImage);
binaryImage = binaryImage+BW*1;
end

binaryImage = binaryImage > 0.5; % Bright objects will be chosen if you use >.
binaryImage = double(binaryImage);
% Size of the image
ImageSizeX = str2double(get(handles.text73,'string'));
ImageSizeY = str2double(get(handles.text78,'string'));

% The ratio of neuron in the image
RatioAllToOne = round(ImageSizeX*ImageSizeY/sum(sum(binaryImage))*10)/10;
set(handles.text81,'string',['1 : ',num2str(RatioAllToOne)]);

figure (88)
subplot(2, 3, 5);
% binaryImage(200:250,1:512)=2;
% binaryImage
imshow(binaryImage); 
title('enlarge the pattern in case of motion', 'FontSize', captionFontSize); 

% binaryImageS = uint8(binaryImage);
% % % 
% binaryImageS = binaryImageS*256;
% imwrite(binaryImageS,'F:\Bo\Resonant_scanner\Imaging_jRGECO_20190329_LongTerm_1_GP858_20190301_M1_\file_00037_Pattern_Structure.tif');


% Convert image to temporal pulses.0
% For 1-direction scanning, such as a polygon  20220323

AWGSamplingRate = str2double(get(handles.edit15,'string'))*1e6;
DTimeAWG = 1/AWGSamplingRate;      % s

ScanBackTime_left     = str2double(get(handles.edit47, 'String'))*DTimeAWG; %Scan back time for left, per line
ScanTime              = str2double(get(handles.edit48, 'String'))*DTimeAWG;
ScanBackTime_Right    = ScanBackTime_left;
ScanBackTime_Right2   = str2double(get(handles.edit50, 'String'))*DTimeAWG;
Scan_Comp_Sample = str2double(get(handles.edit51, 'String'));
FlyBackSamples        = str2double(get(handles.edit52, 'String')); %Flyback time, in samples of AWG.
% ScanAllTime = ScanBackTime_left + ScanTime*2+ScanBackTime_Right+ScanBackTime_Right2;

ScanBackSamples_left  = round(ScanBackTime_left/DTimeAWG); % How many AWG samples for SBL
ScanSamples           = round(ScanTime/DTimeAWG);
ScanBackSamples_Right = round(ScanBackTime_Right/DTimeAWG);    
ScanBackSamples_Right2= round(ScanBackTime_Right2/DTimeAWG);  
ScanAllLineSamples    = ScanBackSamples_left+ScanSamples;
ScanAllLineSamples    = floor(ScanAllLineSamples/4)*4+4;




VTimeSample = linspace(1,ScanAllLineSamples,ScanAllLineSamples);
VTime = DTimeAWG* VTimeSample;

PixelPerLine = ImageSizeX;
NumOfPatterns = ImageSizeY;
PatternLine = zeros(NumOfPatterns-1,ScanAllLineSamples); % 511 for 512;


ImageArea = PatternLine;
IdealRatioOfAllToOnePlusX = str2double(get(handles.edit42,'string'));
UnUsePulseValue = str2double(get(handles.edit39,'string'));
FinalLineRatioAllToOne = zeros(NumOfPatterns,1);
FinalLineRatioAllToOnePulseX = zeros(NumOfPatterns,1);

% One Image: Line 1 to 511
for OI=1:NumOfPatterns-1 % 1:511
    ImageLineFirst  = binaryImage(OI,1:PixelPerLine);
    

    % First scan time, for line 1, 2, 3..., from left to right = 
    % scan time + scan back left (i.e. outside of fill fraction, AWG wait
    % for the next line trigger).
    for FST=1:PixelPerLine
        RPSPixel   = RelationshipPixelAndSampes(FST,1); % e.g., pixel 512
        RPSSamFrom = RelationshipPixelAndSampes(FST,2); % e.g., pixel 512 starts from 1
        RPSSamTo   = RelationshipPixelAndSampes(FST,3); % e.g., pixel 512 ends at 6
        ScanBias   = 0;                                 % number of samples in front of the current image line
        TrailZeros = ScanBackSamples_left;              % number of samples after the current image line
        PatternLine(OI,(RPSSamFrom:RPSSamTo)+ScanBias) = ImageLineFirst(RPSPixel);
        ImageArea(OI,(RPSSamFrom:RPSSamTo)+ScanBias) = 1;
    end
    
    
    % First scan - scanback before it, and gain transient compensation
    SamplesLineAndScanBack = ScanSamples+ScanBackSamples_left;                      % Samples for one line and scan back
    SamplesIdealOneNum = floor(SamplesLineAndScanBack/IdealRatioOfAllToOnePlusX);   % Ideal number of pulses should be on, in one scan.
    ThisScanFrom = RelationshipPixelAndSampes(1,2)+ScanBias;                        % from pixel 512
    ThisScanTo = RelationshipPixelAndSampes(end,3)+ScanBias;                        % to pixel 1
    SamplesActualOneNum = sum(PatternLine(OI,ThisScanFrom:ThisScanTo));             % Actual number of 1 in the first scan
%     if SamplesActualOneNum==0 % New add line, on 20190715 for test,
        if SamplesIdealOneNum>SamplesActualOneNum                                       % if ideal > actual, then need add pulse in right edge, where individual pulses are equaly spaced
            SamplesExtraCompOneNum = SamplesIdealOneNum - SamplesActualOneNum;          % Compensate number of pulses = Ideal - Actual, value will be X, UnUsePulseValue
            RatioOfActualAndComp = SamplesActualOneNum/SamplesExtraCompOneNum;          % Ratio of Actual and compensation number of pulses, i.e., 1:X
            NumOfOneWithKcount  = 0;                                                    % count number of pulses on
            NumOfOneSampleEmployed = 0;                                                 % count number of compensation pulses that are employed
            
%             % comment out below, 20220601 for test
%             for kcount=ThisScanTo:-1:ThisScanFrom                                       % count number of pulse from pixel 1 to pixel 512
%                 NumOfOneWithKcount = NumOfOneWithKcount+PatternLine(OI,kcount);         % count number of pulses
%                 if NumOfOneWithKcount>(RatioOfActualAndComp)                            % if find enough 1 in scan, then add 1 in right scanback
%                     RightkCount = floor(abs(kcount-ThisScanTo)/ScanSamples*ScanBackSamples_Right);% kcount move ScanSamples, then Right kcount move ScanBackSamples_Right
%                     PatternLine(OI,end-RightkCount) = UnUsePulseValue;             % define the compensation pulse value, as UnUsePulseValue
%                     NumOfOneWithKcount = 0;                                             % Once defined the compensation pulse, set count as 0
%                     NumOfOneSampleEmployed = NumOfOneSampleEmployed+1;                  % count the number of compensation pulses employed
%                 end
%             end
%             % comment out above, 20220601 for test
            
            NumOfOneSampleUnemployed = floor(SamplesExtraCompOneNum - NumOfOneSampleEmployed);  % calculate the number of compensation pulses umemployed
            NumOfOneSampleUnemployed = min(NumOfOneSampleUnemployed, TrailZeros);               % if number of one sample unemployed is greater than trailzeros, then just turn on all trailzeros
            VOneSamplePosition = linspace(1,NumOfOneSampleUnemployed,NumOfOneSampleUnemployed); % Define the position of the unemployed pulses
            VOneSamplePosition = round(VOneSamplePosition - NumOfOneSampleUnemployed/2);
            VOneSamplePosition = 0 + VOneSamplePosition + floor(ThisScanTo + TrailZeros/2);        % Put all unemployed pulses in the middle of scanback
            PatternLine(OI,VOneSamplePosition) = UnUsePulseValue;                               % define the unemployed pulse value, as UnUsePulseValue
        else
%             PatternLine(2:(NumOfPatterns-1), floor(ScanBackSamples_left/2))=0.99;    % commented out by SZ 20220110
    %         % Number of 1 should be turned off
    %         SamplesOneShouldBeZero = SamplesOneInRightHalfScan - SamplesIdealOneNum;
    %         RatioOfShould = floor(SamplesOneInRightHalfScan/SamplesOneShouldBeZero);
    %         SOSBNum = 0;
    %         for SOSB=ThisScanFrom:ThisScanTo
    %             if PatternLine(OI,SOSB)==1
    %                 SOSBNum = SOSBNum+1;
    %                 if SOSBNum==RatioOfShould
    % %                         PatternLine(OI,SOSB)=0;
    %                     SOSBNum = 0;
    %                 end
    %             end
    %         end
        end
%     end   %20190715
    
%     PatternLine(1,1:50) = 0.9; 
%     PatternLine(1,60) = 0.99; 
    % Information collect
    SamplesActualOneNum = sum(PatternLine(OI,ThisScanFrom:ThisScanTo) > 0);             % Actual number of 1 in the first scan, re-calculate it
    FinalLineRatioAllToOne(OI)=SamplesLineAndScanBack/SamplesActualOneNum;        % Ratio of (1+x+0):1
    FinalLineRatioAllToOnePulseX(OI)=SamplesLineAndScanBack/SamplesIdealOneNum;   % Ratio of (1+x+0):(1+x)
    
        
    
    
end




% Last: Line 512 and frame flyback
ScanLastLineSamplse = FlyBackSamples+ScanSamples+ScanBackSamples_Right;
ScanLastLineSamplse = floor(ScanLastLineSamplse/4)*4+4;
PatternLastLine     = linspace(0,0,ScanLastLineSamplse);
VTimeSampleLastLine = linspace(1,ScanLastLineSamplse,ScanLastLineSamplse);
ImageAreaLastLine   = PatternLastLine;
for OI=NumOfPatterns
    ImageLineFirst  = binaryImage(OI,1:PixelPerLine);
    % Start with right scanback time

    % First scan time
    for FST=1:PixelPerLine
        RPSPixel   = RelationshipPixelAndSampes(FST,1);% 512
        RPSSamFrom = RelationshipPixelAndSampes(FST,2);%1
        RPSSamTo   = RelationshipPixelAndSampes(FST,3);%6
        ScanBias   = 0;
        TrailZeros = ScanBackSamples_Right;
        PatternLastLine((RPSSamFrom:RPSSamTo)+ScanBias) = ImageLineFirst(RPSPixel);
        ImageAreaLastLine((RPSSamFrom:RPSSamTo)+ScanBias) = 1;
    end

    % Go back to right scanback time
    SamplesLineAndScanBack = ScanSamples+ScanBackSamples_left;                      % Samples for one line and scan back
    SamplesIdealOneNum = floor(SamplesLineAndScanBack/IdealRatioOfAllToOnePlusX);   % Ideal number of pulses should be on, in one scan.
    ThisScanFrom = RelationshipPixelAndSampes(1,2)+ScanBias;                        % from pixel 512
    ThisScanTo = RelationshipPixelAndSampes(PixelPerLine,3)+ScanBias;                        % to pixel 1
    SamplesActualOneNum = sum(PatternLastLine(ThisScanFrom:ThisScanTo));            % Actual number of 1 in the first scan
    if SamplesIdealOneNum>SamplesActualOneNum                                       % if ideal > actual, then need add pulse in right edge
        SamplesExtraCompOneNum = SamplesIdealOneNum - SamplesActualOneNum;          % Compensate number of pulses = Ideal - Actual, value will be X, UnUsePulseValue
        RatioOfActualAndComp = SamplesActualOneNum/SamplesExtraCompOneNum;          % Ratio of Actual and compensation number of pulses, i.e., 1:X
        NumOfOneWithKcount  = 0;                                                    % count number of pulses on
        NumOfOneSampleEmployed = 0;                                                 % count number of compensation pulses that are employed
        for kcount=ThisScanTo:-1:ThisScanFrom                                       % count number of pulse from pixel 1 to pixel 512
            NumOfOneWithKcount = NumOfOneWithKcount+PatternLastLine(kcount);         % count number of pulses
            if NumOfOneWithKcount>(RatioOfActualAndComp)                            % if find enough 1 in scan, then add 1 in right scanback
                RightkCount = floor(abs(kcount-ThisScanTo)/ScanSamples*ScanBackSamples_Right);% kcount move x, then Right kcount move 1
                PatternLastLine(end-RightkCount) = UnUsePulseValue;            % define the compensation pulse value, as UnUsePulseValue
                NumOfOneWithKcount = 0;                                             % Once defined the compensation pulse, set count as 0
                NumOfOneSampleEmployed = NumOfOneSampleEmployed+1;                  % count the number of compensation pulses employed
            end
        end
        NumOfOneSampleUnemployed = floor(SamplesExtraCompOneNum - NumOfOneSampleEmployed);  % calculate the number of compensation pulses umemployed
        NumOfOneSampleUnemployed = min(NumOfOneSampleUnemployed, TrailZeros);
        VOneSamplePosition = linspace(1,NumOfOneSampleUnemployed,NumOfOneSampleUnemployed); % Define the position of the unemployed pulses
        VOneSamplePosition = round(VOneSamplePosition - NumOfOneSampleUnemployed/2);
        VOneSamplePosition = 0 + VOneSamplePosition+ floor(ThisScanTo + TrailZeros/2);                                        % Define the position of the unemployed pulses
        PatternLastLine(VOneSamplePosition) = UnUsePulseValue;                               % define the unemployed pulse value, as UnUsePulseValue
        
       
    else
%         % Number of 1 should be turned off
%         SamplesOneShouldBeZero = SamplesOneInRightHalfScan - SamplesIdealOneNum;
%         RatioOfShould = floor(SamplesOneInRightHalfScan/SamplesOneShouldBeZero);
%         SOSBNum = 0;
%         for SOSB=ThisScanFrom:ThisScanTo
%             if PatternLine(OI,SOSB)==1
%                 SOSBNum = SOSBNum+1;
%                 if SOSBNum==RatioOfShould
% %                         PatternLine(OI,SOSB)=0;
%                     SOSBNum = 0;
%                 end
%             end
%         end
    end
%      PatternLastLine(floor(ScanBackSamples_left/2))=0.99;
    % Information collect
    SamplesActualOneNum = sum(PatternLastLine(ThisScanFrom:ThisScanTo) > 0); % Actual number of 1 in half center scan
    FinalLineRatioAllToOne(end)=SamplesLineAndScanBack/SamplesActualOneNum;
    FinalLineRatioAllToOnePulseX(end)=SamplesLineAndScanBack/SamplesActualOneNum;

    % Left scan back time, 4 MHz
    FlyBackSamFrom = RPSSamTo+ScanBias+8;
    FlyBackSamTo = RPSSamTo+ScanBias+FlyBackSamples-1;
%     disp('pattern last line length 1: ' + string(length(PatternLastLine)));
    PatternLastLine(round(FlyBackSamFrom:round(IdealRatioOfAllToOnePlusX/str2double(get(handles.edit41,'string'))):FlyBackSamTo)) = UnUsePulseValue;

    
    

    % Second scan time
%     for SST=1:512
%         RPSPixel   = RelationshipPixelAndSampes(SST,1);% 512
%         RPSSamFrom = RelationshipPixelAndSampes(SST,2);%1
%         RPSSamTo   = RelationshipPixelAndSampes(SST,3);%6
%         ScanBias = ScanBackSamples_Right+ScanSamples+FlyBackSamples;
%         PatternLastLine((RPSSamFrom:RPSSamTo)+ScanBias) = ImageLineSecond(513-RPSPixel);
%         ImageAreaLastLine((RPSSamFrom:RPSSamTo)+ScanBias) = 1;
%     end

    % One line finished

end
PatternLine(1,2:3)=0;
figure (88)
subplot(4, 3, 4);
for k=1
    VTimeTem = VTimeSample+ScanAllLineSamples*(k-1);
    plot(VTimeTem,PatternLine(k,:),'b');hold on
    plot(VTimeTem,ImageArea(k,:),'-.r');hold on
end
xlim([min(VTimeTem),max(VTimeTem)]);
ylim([-0.05,1.05]);
title('Pattern in time, Line 1', 'FontSize', captionFontSize); 
% xlabel('Time/AWG Samples');

figure (88)
subplot(4, 3, 7);
for k=2
    VTimeTem = VTimeSample+ScanAllLineSamples*(k-1);
    plot(VTimeTem,PatternLine(k,:),'b');hold on
    plot(VTimeTem,ImageArea(k,:),'-.r');hold on
end
xlim([min(VTimeTem),max(VTimeTem)]);
ylim([-0.05,1.05]);
% xlabel('Time/AWG Samples');
title('Pattern in time, Line 2', 'FontSize', captionFontSize); 
handles.ConvertedPatternInSamples = PatternLine;
% PatternLineLastOne = 
figure (88)
subplot(4, 3, 10);

    VTimeTem = VTimeSampleLastLine+ScanAllLineSamples*(NumOfPatterns-1);
    length(PatternLastLine)
    length(VTimeTem)
    plot(VTimeTem,PatternLastLine,'b');hold on
    plot(VTimeTem,ImageAreaLastLine,'-.r');hold on

xlim([min(VTimeTem),max(VTimeTem)]);
ylim([-0.05,1.05]);
xlabel('Time/AWG Samples');
title('Pattern in time, last line', 'FontSize', captionFontSize); 
handles.ConvertedPatternInSamplesFirst = PatternLine;
% csvwrite('firstPattern.csv', PatternLine);
whos PatternLine
handles.ConvertedImageArea = ImageArea;
handles.ConvertedPatternInSamplesLast = PatternLastLine;
% csvwrite('lastPattern.csv', PatternLastLine);


    
    
    % First part, line 2 to 511
    PatternInSampleFirst = PatternLine;
    LineOfThePatternFirst = length(PatternInSampleFirst(:,1));
    SamplesOfOneLineFirst = length(PatternInSampleFirst(1,:));
    SamplesOfOnePatternFirst = (SamplesOfOneLineFirst+Scan_Comp_Sample)*LineOfThePatternFirst; % Add Scan_Comp_Sample samples to two lines: for example, line 2 and 3.
    NumOfOneInFirst = sum(sum(PatternInSampleFirst>0));
    RatioOfZeroToOne = round(SamplesOfOnePatternFirst/NumOfOneInFirst*10)/10;
    % Second part, line 512 to line 1
    PatternInSampleLast = PatternLastLine;
    SamplesOfOneLineLast = length(PatternInSampleLast);
    
    SamplesOfOnePattern = SamplesOfOnePatternFirst + SamplesOfOneLineLast;
    SamplesOfOnePattern=floor(SamplesOfOnePattern/4)*4+4;
    NSampleAWG = SamplesOfOnePattern;
    TWindowAWG = DTimeAWG*NSampleAWG;
    
    set(handles.text88,'string',['1:',num2str(RatioOfZeroToOne)]);
    set(handles.text97,'string',num2str(TWindowAWG/1e-6));
    set(handles.edit16,'string',num2str(floor(TWindowAWG/1e-6+100)));
    set(handles.text94,'string','Yes');

% Image information
PixelSumPerline = sum(binaryImage');
PixelSumPerlineMax = max(PixelSumPerline);
set(handles.text100,'string',[num2str(PixelSumPerlineMax),',1:',num2str(round(PixelPerLine/PixelSumPerlineMax*10)/10)]);
% set(handles.text102,'string',[num2str(PixelSumPerlineMin)]);

% Pattern information
SamplesSumPerline = sum((PatternLine(2:end,:)>0),2);
SamplesSumPerlineMax = max(SamplesSumPerline);

% 
% PixelSumPerlineMin = min(PixelSumPerline);
set(handles.text102,'string',['1:',num2str(round(ScanAllLineSamples/SamplesSumPerlineMax*10)/10)]);
% set(handles.text102,'string',[num2str(PixelSumPerlineMin)]);


% Information of image

BigNumber = 10000;
FinalLineRatioAllToOne(end)=BigNumber;
FinalLineRatioAllToOnePulseX(end)=BigNumber;

minFinalLineRatioAllToOne = min(FinalLineRatioAllToOne);
AveFinalLineRatioAllToOne = (sum(1./FinalLineRatioAllToOne)-1/BigNumber)/(NumOfPatterns);

minFinalLineRatioAllToOnePulseX = min(FinalLineRatioAllToOnePulseX);
AveFinalLineRatioAllToOnePulseX = (sum(1./FinalLineRatioAllToOnePulseX)-1/BigNumber)/(NumOfPatterns);
% whos FinalLineRatioAllToOnePulseX

set(handles.text48 ,'string',['max = 1:',num2str(round(minFinalLineRatioAllToOne*10)/10),', average = 1:',num2str(round(1/AveFinalLineRatioAllToOne*10)/10)]);
set(handles.text129,'string',['max = 1:',num2str(round(minFinalLineRatioAllToOnePulseX*10)/10),', average = 1:',num2str(round(1/AveFinalLineRatioAllToOnePulseX*10)/10)]);


        

% VIntensityAWG         = linspace(PulseOff,PulseOff,NSampleAWG); % Pre-dine AWG pattern
% VIntensityAWG(560+1) = PulseOn;
% VIntensityAWG(600:1000) = PulseOn;
% VIntensityAWG(2600:2700) = PulseOn;
% 
% 
% 
% 
% LineToAllScanRatio = str2double(get(handles.edit39,'string'));
% PixelPerLine = ImageSizeX;
% PixelPerScan = PixelPerLine/LineToAllScanRatio;
% PixelUnusePerLine = PixelPerScan - PixelPerLine;
% LinePerFrame = ImageSizeY;
% 
% % calculate the sum of pixel values in every line, and the maximum sum
% PixelsEveryLine = linspace(0,0,LinePerFrame);
% for kpixelperline=1:LinePerFrame
%     PixelsEveryLine(kpixelperline) = sum(binaryImage(kpixelperline,:));
%     MaxPixelOnPerLine = max(PixelsEveryLine);
%     MinPixelOnPerLine = min(PixelsEveryLine);
% end
% set(handles.text100,'string',num2str(MaxPixelOnPerLine));
% set(handles.text102,'string',num2str(MinPixelOnPerLine));
% 
% % Set the pixel time
% PixelTime = str2double(get(handles.edit42,'string'))*1e-6;
% 
% % number of pixels
% NTime = floor(PixelPerLine*LinePerFrame/LineToAllScanRatio)+1;
% % pixel vector
% VPixels = linspace(1,NTime,NTime);
% % Time vector
% VTime = VPixels*PixelTime;
% % Time vector for scanning a line (include useless area)
% TimePerScan = PixelPerScan * PixelTime;
% % Time vector for scanning a line (only consists of useful image area)
% TimePerLine = PixelPerLine * PixelTime;
% % Intensity of pixel vector
% VIntensity = VTime*0+0;
% 
% VLineArea = VTime*0;
% CompForBlankLine = str2double(get(handles.edit40,'string')); % to avoid no pulse in a line.
% CompStartRatio   = str2double(get(handles.edit41,'string'))/100; % where does the compensation starts.
% % it is the ratio of start to all pixels that are not used.
% EqualEachLine = get(handles.text104,'string'); % to avoid no pulse in a line.
% % binaryImage
% % if yes, then equal each line. if no, just make sure every line is not
% % totally off.
% if EqualEachLine == 'Yes'
%     for k=1:LinePerFrame
%         StartPoint = round(1+PixelPerScan*(k-1));
%         VIntensity(StartPoint:StartPoint+PixelPerLine-1) = binaryImage(k,1:PixelPerLine);
%         VLineArea(StartPoint:StartPoint+PixelPerLine-1) = 1;
%         % Making sure every line has the same on pulses.
% %         StarPoint2 = StartPoint+PixelPerLine+floor(PixelUnusePerLine*CompStartRatio);
% %         VIntensity(StarPoint2:StarPoint2+MaxPixelOnPerLine-PixelsEveryLine(k)-1) = 1;    
% %         StarPoint2 = StartPoint+PixelPerLine+floor(PixelUnusePerLine*CompStartRatio)*0;
% %                 CompensationNumber = floor((MaxPixelOnPerLine-PixelsEveryLine(k))*0.5);
% %         VIntensity(StarPoint2:StarPoint2+CompensationNumber-1) = 0.9;    
% %         VIntensity(StartPoint+PixelPerScan-CompensationNumber:StartPoint+PixelPerScan-1) = 0.9;    
%         StarPoint2 = StartPoint+PixelPerLine+floor(PixelUnusePerLine*CompStartRatio);
%         VIntensity(StarPoint2:StarPoint2+MaxPixelOnPerLine-PixelsEveryLine(k)-1) = 1;    
%     end
% else
%     for k=1:LinePerFrame
%         StartPoint = round(1+PixelPerScan*(k-1));
%         VIntensity(StartPoint:StartPoint+PixelPerLine-1) = binaryImage(k,1:PixelPerLine);
%         VLineArea(StartPoint:StartPoint+PixelPerLine-1) = 1;
%         % Extra pulses to ensure some pulses in every line
%         if sum(binaryImage(k,:))<=1
%             StarPoint2 = StartPoint+PixelPerLine+floor(PixelUnusePerLine*CompStartRatio);
%             VIntensity(StarPoint2:StarPoint2+CompForBlankLine-1) = 1;
%         end
%         % add pulse in the beginning
% %         
% %             StarPoint2 = StartPoint+PixelPerLine+floor(PixelUnusePerLine*CompStartRatio);
% %             VIntensity(StarPoint2:StarPoint2+CompForBlankLine-1) = 1;
%     end    
% end
% 
% 
% 
% % Final pulse on to all, ratio
% FinalPulseAlltoOnRatio = round( NTime/sum(VIntensity)*10 )/10;
% set(handles.text88,'string',['1 : ',num2str(FinalPulseAlltoOnRatio)]);
% 
% AWGPatternExtraLines = str2double(get(handles.edit44,'string'));
% AWGPatternExtraPixels = AWGPatternExtraLines*PixelPerScan;
% AWGPatternExtraPixelsOn = AWGPatternExtraPixels.*FinalPulseAlltoOnRatio;
% 
% 
% 
% figure (88)
% subplot(4, 3, 4);
% plot(VTime,VLineArea,'-.r');hold on
% hold on;
% % whos VTime
% % whos VIntensity
% plot(VTime,VIntensity,'b');
% 
% hold on;
% ylim([-0.2,1.6]);
% xlim([min(VTime),max(VTime)]);
% ylabel('Intensity/a.u.');
% title('Converted pattern in time domain, in second', 'FontSize', captionFontSize); 
% 
% subplot(4, 3, 7);
% plot(VTime/1e-6,VIntensity,'b');
% hold on;
% plot(VTime/1e-6,VLineArea,'-.r');
% hold on;
% ylim([-0.2,1.8]);
% xlim([0,TimePerScan/1e-6*2]);
% xlabel('Time/us');
% ylabel('Intensity/a.u.');
% % Line 1
% text(0,1.3,'------Line 1-------','Color','r');
% text(TimePerScan/1e-6,1.3,'------Line 2-------','Color','r');
% text(TimePerLine/1e-6,1.6,'-edge of scan-','Color','k');
% text(TimePerLine/1e-6+TimePerScan/1e-6,1.6,'-edge of scan-','Color','k');
% 
% subplot(4, 3, 10);
% plot(VPixels,VIntensity,'b');
% hold on;
% plot(VPixels,VLineArea,'-.r');
% hold on;
% ylim([-0.2,1.8]);
% xlim([0,PixelPerScan*2]);
% xlabel('Pixels');
% ylabel('Intensity/a.u.');
% % Line 1
% text(0,1.3,'------Line 1-------','Color','r');
% text(PixelPerScan,1.3,'------Line 2-------','Color','r');
% text(PixelPerLine,1.6,'-edge of scan-','Color','k');
% text(PixelPerLine+PixelPerScan,1.6,'-edge of scan-','Color','k');
% 
% handles.ConvertedPatternInPixel = VIntensity;
% set(handles.text94,'string','Yes');
% 
% 
% AWGSamplingRate = str2double(get(handles.edit15,'string'))*1e6;
% DTimeAWG = 1/AWGSamplingRate;      % s
% AWGSamplesPerPixel = (PixelTime/DTimeAWG);
% NPatternInPixel = length(VIntensity);
% NSampleAWG = floor(NPatternInPixel*AWGSamplesPerPixel); % AWG has equal or more samples than pixels.
% NSampleAWG = floor(NSampleAWG/4)*4+4;     % Sample number should be 4X
% TWindowAWG = DTimeAWG*NSampleAWG;
% set(handles.text97,'string',num2str(TWindowAWG/1e-6));
% set(handles.edit16,'string',num2str(floor(TWindowAWG/1e-6+100)));

Newhandles = handles;

% figure
% plot(VTime,VLineArea,'-.r');
% hold on;
% plot(VTime,VIntensity,'k');
% hold on;
% ylim([-0.5,1.5]);



% tt=linspace(1,512,512);
% % binaryImage = binaryImage/255
% Line1 = binaryImage(260,:);
% Line1Pro =double(Line1'>0);
% Ratio1 = sum(Line1Pro)/512*100;
% Line2 = binaryImage(175,:);
% Line2Pro =double(Line2'>0);
% Ratio2 = sum(Line2Pro)/512*100;
% Line3 = binaryImage(395,:);
% Line3Pro =double(Line3'>0);
% Ratio3 = sum(Line3Pro)/512*100;
% LineAll = [Line1Pro,Line2Pro,Line3Pro];
% figure (88)
% subplot(2, 3, 4);
% plot(tt,Line2Pro,'linewidth',3);hold on;
% plot(tt,Line1Pro+1.5,'linewidth',3);hold on;
% plot(tt,Line3Pro+3,'linewidth',3);hold on;
% ylim([-0.5,4.5]);
% xlim([1,512]);
% xlabel('Time/ *pixel time');
% ylabel('Intensity/a.u.')

% 
% % % Display the processed binary image.
% % binaryImage = bwareaopen(binaryImage,60);
% % subplot(3, 3, 5);
% % imshow(binaryImage); 
% % title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize); 
% 
% % Identify individual blobs by seeing which pixels are connected to each other.
% % Each group of connected pixels will be given a label, a number, to identify it and distinguish it from the other blobs.
% % Do connected components labeling with either bwlabel() or bwconncomp().
% labeledImage = bwlabel(binaryImage, 8);     % Label each blob so we can make measurements of it
% % labeledImage is an integer-valued image where all pixels in the blobs have values of 1, or 2, or 3, or ... etc.
% figure (88)
% subplot(3, 3, 8);
% imshow(labeledImage, []);  % Show the gray scale image.
% % title('Labeled Image, from bwlabel()', 'FontSize', captionFontSize);
% 
% 
% 
% % Let's assign each blob a different color to visually show the user the distinct blobs.
% coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
% % coloredLabels is an RGB image.  We could have applied a colormap instead (but only with R2014b and later)
% figure (88)
% subplot(3, 3, 9);
% imshow(coloredLabels);
% axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
% % caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
% title(caption, 'FontSize', captionFontSize);
% 
% % Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
% blobMeasurements = regionprops(labeledImage, originalImage, 'all');
% numberOfBlobs = size(blobMeasurements, 1);
% 
% 
% % bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
% % Plot the borders of all the coins on the original grayscale image using the coordinates returned by bwboundaries.
% figure (88)
% subplot(3, 3, 7);
% imshow(originalImage);
% title('Outlines', 'FontSize', captionFontSize); 
% axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
% hold on;
% boundaries = bwboundaries(binaryImage);
% numberOfBoundaries = size(boundaries, 1);
% for k = 1 : numberOfBoundaries
% 	thisBoundary = boundaries{k};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
% end
% hold off;
% 
% textFontSize = 14;	% Used to control size of "blob number" labels put atop the image.
% labelShiftX = -7;	% Used to align the labels in the centers of the coins.
% blobECD = zeros(1, numberOfBlobs);
% % Print header line in the command window.
% % fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
% % Loop over all blobs printing their measurements to the command window.
% for k = 1 : numberOfBlobs           % Loop through all blobs.
% 	% Find the mean of each blob.  (R2008a has a better way where you can pass the original image
% 	% directly into regionprops.  The way below works for all versions including earlier versions.)
% 	thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
% 	meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
% 	meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
% 	
% 	blobArea = blobMeasurements(k).Area;		% Get area.
% 	blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
% 	blobCentroid = blobMeasurements(k).Centroid;		% Get centroid one at a time
% 	blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
% % 	fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));
% 	% Put the "blob number" labels on the "boundaries" grayscale image.
% 	text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold');
% end

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImageProcessed = ImageProcessStackProcess(hObject, eventdata, handles);
% Maximize the figure window.
figure
imshow(ImageProcessed);
caption = sprintf('Fluctuation of each pixel in the stack');
title(caption, 'FontSize', handles.captionFontSize);
axis image; % Make sure image is not artificially stretched because of screen's

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fighandles = findall( allchild(0), 'type', 'figure');
NumOfFigs = length(fighandles);
for k=1:NumOfFigs
    Number = fighandles(k).Number;
    if isempty(Number)~=1
    delete(Number);
    end
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text44, 'String', 'Channel A'); 

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text44, 'String', 'Channel B'); 

% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text44, 'String', 'External '); 



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AWG8on8off(hObject, eventdata, handles);


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% PresentText = get(handles.text48, 'String'); 
% if PresentText == 'Yes'
%     set(handles.text48, 'String', 'No ');
% else
%     set(handles.text48, 'String', 'Yes');
% end


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text50, 'String', 'A  ');

% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text50, 'String', 'B  ');

% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text50, 'String', 'A+B');



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles]=OSCMeasureAndSave(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);
% OSC Initialization status update
OSCInitializationCheck(handles);






function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentText = get(handles.text70, 'String'); 
if PresentText == 'Stack'
    set(handles.text70, 'String', 'Image');
else
    set(handles.text70, 'String', 'Stack');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentText = get(handles.text94, 'String'); 
if PresentText == 'Yes'
    set(handles.text94, 'String', 'No ');
else
    set(handles.text94, 'String', 'Yes');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentText = get(handles.text104, 'String'); 
if PresentText == 'Yes'
    set(handles.text104, 'String', 'No ');
else
    set(handles.text104, 'String', 'Yes');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentText = get(handles.text108, 'String'); 
if PresentText == 'Yes'
    set(handles.text108, 'String', 'No ');
else
    set(handles.text108, 'String', 'Yes');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AWG32MHzTo4MHz(hObject, eventdata, handles);



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StepOrContinuous = 1; % Stepped
RunLastWaveformStepped(hObject, handles, StepOrContinuous);


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[hObject, handles]=AWGDisable(hObject, handles);


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StepOrContinuous = 0; % Continuous
RunLastWaveform(hObject, handles, StepOrContinuous);



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

StepOrContinuous = 1; % Stepped
RunLastWaveform(hObject, handles, StepOrContinuous);



function edit53_Callback(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit53 as text
%        str2double(get(hObject,'String')) returns contents of edit53 as a double


% --- Executes during object creation, after setting all properties.
function edit53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PresentText = get(handles.text134, 'String'); 
if PresentText == 'Internal'
    set(handles.text134, 'String', 'External');
    set(handles.edit11, 'String', num2str(672));
    set(handles.text50, 'String', 'A  ');    
else
    set(handles.text134, 'String', 'Internal');
    set(handles.edit11, 'String', num2str(1000));
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit56_Callback(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit56 as text
%        str2double(get(hObject,'String')) returns contents of edit56 as a double


% --- Executes during object creation, after setting all properties.
function edit56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit57_Callback(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit57 as text
%        str2double(get(hObject,'String')) returns contents of edit57 as a double


% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit58_Callback(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58 as text
%        str2double(get(hObject,'String')) returns contents of edit58 as a double


% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.text144, 'String', 'Simple');

% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text144, 'String', 'Full  ');

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text144, 'String', 'NO    ');
