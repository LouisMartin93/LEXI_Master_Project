% this is a test script fo initialisation

clear; close all; clc;
global h1; % make h a global variable so it can be used outside the main
global h2; % function. Useful when you do event handling and sequential           move
%% Create Matlab Figure Container
fpos    = get(0,'DefaultFigurePosition'); % figure default position
fpos(1) = 10; % figure window size;Width
fpos(2) = 290;
fpos(3) = 650; % figure window size;Width
fpos(4) = 450; % Height

f1 = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI');
       

fpos(1) = 675; % figure window size;Width
fpos(2) = 290; % Height
f2 = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI'); 

%% Create ActiveX Controller
h1 = actxcontrol('MGMOTOR.MGMotorCtrl.1',[20 20 600 400 ], f1);
h1.StartCtrl;
% Set the Serial Number
SN = 83817744; % put in the serial number of the hardware
set(h1,'HWSerialNum', SN);
h1.Identify;
pause(5);

%% Create ActiveX Controller
h2 = actxcontrol('MGMOTOR.MGMotorCtrl.1',[20 20 600 400 ], f2);
h2.StartCtrl;
% Set the Serial Number
SN = 83817899; % put in the serial number of the hardware
set(h2,'HWSerialNum', SN);
h2.Identify;
pause(5);
%% home position both of these
h1.MoveHome(0,0);
h2.MoveHome(0,0);


