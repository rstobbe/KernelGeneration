%=========================================================
% 
%=========================================================

function [KRNprms,err] = Design_ConvKernel_v1a_Func(INPUT,KRNprms)

Status('busy','Create Kernel');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Get Input
%---------------------------------------------
KDES = INPUT.KDES;
clear INPUT;

%---------------------------------------------
% Run
%---------------------------------------------
func = str2func([KRNprms.kerndesfunc,'_Func']);
INPUT = [];
[KDES,err] = func(KDES,INPUT);
if err.flag
    return
end
clear INPUT

%---------------------------------------------
% Add to Panel Output
%---------------------------------------------
Panel(1,:) = {'','','Output'};
Panel(2,:) = {'',KRNprms.method,'Output'};

KRNprms = KDES;
KRNprms.Panel = [Panel;KRNprms.Panel];
KRNprms.PanelOutput = cell2struct(KRNprms.Panel,{'label','value','type'},2);
KRNprms.ExpDisp = PanelStruct2Text(KRNprms.PanelOutput);


Status('done','');
Status2('done','',2);
Status2('done','',3);
