% Set up script that will run to set up the camera make a master dark

clear; close all;

%adding the external paths for the various scripts and functions
addpath('Motor_ctrl');
addpath('PIKE_Matlab');
addpath('spot extraction noGUI');


global h1; % make h a global variable so it can be used outside the main
global h2; % function. Useful when you do event handling and sequential moves

%% Create Matlab Figure Container
fpos    = get(0,'DefaultFigurePosition'); % figure default position
fpos(1) = 10; % figure window size;Width
fpos(2) = 290;
fpos(3) = 650; % figure window size;Width
fpos(4) = 450; % Height

f1 = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI 1');
       

fpos(1) = 675; % figure window size;Width
fpos(2) = 290; % Height
f2 = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI 2'); 

%% Create ActiveX Controller 1
SN1 = 83817744;
h1 = CtrlSetup('MGMOTOR.MGMotorCtrl.1',SN1,f1);

%% Create ActiveX Controller 2
SN2 = 83817899;
h2 = CtrlSetup('MGMOTOR.MGMotorCtrl.1',SN2,f2);

%% Setting to Motors to 'home' 
% This will move these to home. The calculated 180 and 0 position of P1 and
% P2 are calculated as 153.4 and 6.01
ADCMovePos(h1,153.4);
ADCMovePos(h2,6.01);
% now the prisms are at the canceling position so that the beam passes
% straight through.


%% Initialising the camera can be one manually with user input or automated

%manually
vid = init_camera;

%{
%automatic
vid = videoinput('dcam', 1,'F7_Y16_640x480');
vid.FramesPerTrigger = 1;
src = getselectedsource(video_object);
src.Shutter = 1000; %(us)
%}

%% Creating the Dark image
DarkDir = 'PIKE_Matlab\Test Images\dispersion\Setup 1\SM450\f4_cMWS\darks\';
MasterDark = makeMasterDark(DarkDir);
pause(41.4) %pause for motor setup
disp('Setup Complete')