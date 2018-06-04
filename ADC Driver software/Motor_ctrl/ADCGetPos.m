function position = ADCGetPos(Ctrl_Object)
    position = Ctrl_Object.GetPosition_Position(0);
    disp(position);
end