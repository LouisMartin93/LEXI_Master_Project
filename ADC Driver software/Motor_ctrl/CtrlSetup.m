function [Ctrl_Object]=CtrlSetup(ActXControlID,SerialNumber,fig)
    % defining the contorl object
    Ctrl_Object = actxcontrol(ActXControlID,[20 20 600 400],fig);
    % starting communication between computer and object
    Ctrl_Object.StartCtrl;
    % set the serial number of the hardware
    set(Ctrl_Object,'HWSerialNum',SerialNumber);
    % Identify the Object
    Ctrl_Object.Identify;
    pause(5); % to allow identification
    addproperty(Ctrl_Object,'FIGURE');
    Ctrl_Object.set('FIGURE',fig.Number);
    %Moves the motor to its home position
    Ctrl_Object.MoveHome(0,0);
end
    