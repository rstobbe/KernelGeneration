%==================================================
% (v1d) 
%       - User ParFor
%       - Convscaleval resolution update
%==================================================

function [SCRPTipt,SCRPTGBL,err] = SglC_KB3DCart_v1d(SCRPTipt,SCRPTGBL)

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
KRNprms.W = str2double(SCRPTGBL.CurrentScript.('Width'));
KRNprms.res = str2double(SCRPTGBL.CurrentScript.('Res'));
KRNprms.beta = str2double(SCRPTGBL.CurrentScript.('Beta'));
KRNprms.DesforSS = str2double(SCRPTGBL.CurrentScript.('DesForSS')); 
KRNprms.prec = 'singles';

%---------------------------------------------
% Create Kernel
%---------------------------------------------
func = str2func([KRNprms.method,'_Func']);
INPUT = [];
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
SCRPTipt(indnum).entrystr = cell2mat(name);

SCRPTGBL.RWSUI.SaveVariables = {KRNprms};
SCRPTGBL.RWSUI.SaveVariableNames = 'KRNprms';
SCRPTGBL.RWSUI.SaveGlobal = 'yes';
SCRPTGBL.RWSUI.SaveGlobalNames = name;
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = 'outloc';
SCRPTGBL.RWSUI.SaveScriptName = cell2mat(name);

Status('done','');
Status2('done','',2);
Status2('done','',3);

