%==================================================
% (v1b)
%       - as Inv_KB3DCart_v1b (save as single)
%==================================================

function [SCRPTipt,SCRPTGBL,err] = Inv_KB3DCartS_v1b(SCRPTipt,SCRPTGBL)

Status('busy','Generate Inverse Filter');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Clear Naming
%---------------------------------------------
inds = strcmp('IF_Name',{SCRPTipt.labelstr});
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
IF.script = SCRPTGBL.CurrentTree.Func;
IF.W = str2double(SCRPTGBL.CurrentScript.Width);
IF.beta = str2double(SCRPTGBL.CurrentScript.Beta);
IF.SS = str2double(SCRPTGBL.CurrentScript.UseWithSS); 
IF.ZF = str2double(SCRPTGBL.CurrentScript.UseWithZF);
IF.anlzres = str2double(SCRPTGBL.CurrentScript.anlzRes);
IF.anlzzf = str2double(SCRPTGBL.CurrentScript.anlzZF);

%---------------------------------------------
% Generate
%---------------------------------------------
func = str2func([IF.script,'_Func']);
INPUT.IF = IF;
[OUTPUT,err] = func(INPUT);
if err.flag
    return
end

%--------------------------------------------
% Output
%--------------------------------------------
IFprms = OUTPUT.IFprms;
clear OUTPUT

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
name = inputdlg('Name IF Filter:','Name IF Filter',1,{IFprms.name});
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

