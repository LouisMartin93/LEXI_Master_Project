function ADCMovePos(Ctrl_Object,position)
    Ctrl_Object.SetAbsMovePos(0,position);
    Ctrl_Object.MoveAbsolute(0,0==1);
end