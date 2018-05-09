%set up scrit v2, change is that set up now uses functions

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

