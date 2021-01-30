%==================================================
% (v1c)
%       - Interface Fixup
%==================================================

function [SCRPTipt,SCRPTGBL,err] = Inv_KB3DCartS_v1c(SCRPTipt,SCRPTGBL)

Status('busy','Generate Inverse Filter');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Clear Naming
%---------------------------------------------
inds = strcmp('InvFilt_Name',{SCRPTipt.labelstr});
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
IFprms.method = SCRPTGBL.CurrentTree.Func;
IFprms.W = str2double(SCRPTGBL.CurrentScript.Width);
IFprms.beta = str2double(SCRPTGBL.CurrentScript.Beta);
IFprms.SS = str2double(SCRPTGBL.CurrentScript.UseWithSS); 
IFprms.ZF = str2double(SCRPTGBL.CurrentScript.UseWithZF);
IFprms.anlzres = str2double(SCRPTGBL.CurrentScript.anlzRes);
IFprms.anlzzf = str2double(SCRPTGBL.CurrentScript.anlzZF);

%---------------------------------------------
% Generate
%---------------------------------------------
func = str2func([IFprms.method,'_Func']);
INPUT = struct();
[IFprms,err] = func(INPUT,IFprms);
if err.flag
    return
end

%--------------------------------------------
% Output to TextBox
%--------------------------------------------
% DES.ExpDisp = PanelStruct2Text(DES.PanelOutput);
% DES.ExpDisp = [newline DES.ExpDisp];
% global FIGOBJS
% FIGOBJS.(SCRPTGBL.RWSUI.tab).Info.String = DES.ExpDisp;

%--------------------------------------------
% Return
%--------------------------------------------
name = inputdlg('Name Inverse Filter:','Name Inverse Filter',1,{IFprms.name});
if isempty(name)
    SCRPTGBL.RWSUI.SaveGlobal = 'no';
    return
end
IFprms.name = name{1};

SCRPTipt(indnum).entrystr = IFprms.name;
SCRPTGBL.RWSUI.SaveVariables = IFprms;
SCRPTGBL.RWSUI.SaveVariableNames = 'IFprms';
SCRPTGBL.RWSUI.SaveGlobal = 'no';
SCRPTGBL.RWSUI.SaveGlobalNames = IFprms.name;
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = 'outloc';
SCRPTGBL.RWSUI.SaveScriptName = IFprms.name;

Status('done','');
Status2('done','',2);
Status2('done','',3);

