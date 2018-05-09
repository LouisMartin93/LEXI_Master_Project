function CtrlShutdown(Ctrl_Object)
    FigNum = Ctrl_Object.get('FIGURE');
    Ctrl_Object.StopCtrl;
    close(FigNum);
    clear Ctrl_Object
end