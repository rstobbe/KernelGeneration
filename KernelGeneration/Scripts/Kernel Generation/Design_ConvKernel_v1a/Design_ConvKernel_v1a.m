%==================================================
% (v1a) 
%       
%==================================================

function [SCRPTipt,SCRPTGBL,err] = Design_ConvKernel_v1a(SCRPTipt,SCRPTGBL)

Status('busy','Create Convolution Kernel');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Clear Naming
%---------------------------------------------
inds = strcmp('Kernel_Name',{SCRPTipt.labelstr});
indnum = find(inds==1);
if length(indnum) > 1
    indnum = indnum(SCRPTGBL.RWSUI.scrptnum);
end
SCRPTipt(indnum).entrystr = '';
setfunc = 1;
DispScriptParam(SCRPTipt,setfunc,SCRPTGBL.RWSUI.tab,SCRPTGBL.RWSUI.panelnum);

%---------------------------------------------
% Load Input
%---------------------------------------------
KRNprms.method = SCRPTGBL.CurrentTree.Func;
KRNprms.kerndesfunc = SCRPTGBL.CurrentTree.('Kernfunc').Func;

%---------------------------------------------
% Get Working Structures from Sub Functions
%---------------------------------------------
KDESipt = SCRPTGBL.CurrentTree.('Kernfunc');
if isfield(SCRPTGBL,('Kernfunc_Data'))
    KDESipt.Kernfunc_Data = SCRPTGBL.Kernfunc_Data;
end

%------------------------------------------
% Get Function Info
%------------------------------------------
func = str2func(KRNprms.kerndesfunc);           
[SCRPTipt,KDES,err] = func(SCRPTipt,KDESipt);
if err.flag
    return
end

%---------------------------------------------
% Create Kernel
%---------------------------------------------
func = str2func([KRNprms.method,'_Func']);
INPUT.KDES = KDES;
[KRNprms,err] = func(INPUT,KRNprms);
if err.flag
    return
end

%--------------------------------------------
% Output to TextBox
%--------------------------------------------
KRNprms.ExpDisp = PanelStruct2Text(KRNprms.PanelOutput);
global FIGOBJS
FIGOBJS.(SCRPTGBL.RWSUI.tab).Info.String = KRNprms.ExpDisp;

%--------------------------------------------
% Return
%--------------------------------------------
name = inputdlg('Name Kernel:','Name',1,{KRNprms.name});
if isempty(name)
    SCRPTGBL.RWSUI.SaveGlobal = 'no';
    return
end
KRNprms.name = name{1};

SCRPTipt(indnum).entrystr = KRNprms.name;
SCRPTGBL.RWSUI.SaveVariables = KRNprms;
SCRPTGBL.RWSUI.SaveVariableNames = 'KRNprms';
SCRPTGBL.RWSUI.SaveGlobal = 'yes';
SCRPTGBL.RWSUI.SaveGlobalNames = KRNprms.name;
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = 'outloc';
SCRPTGBL.RWSUI.SaveScriptName = KRNprms.name;

Status('done','');
Status2('done','',2);
Status2('done','',3);

