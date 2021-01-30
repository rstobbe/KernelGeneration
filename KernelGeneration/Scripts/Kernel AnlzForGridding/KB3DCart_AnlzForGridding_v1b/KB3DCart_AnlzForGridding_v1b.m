%==================================================
% (v1a)
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = KB3DCart_AnlzForGridding_v1b(SCRPTipt,SCRPTGBL)

Status('busy','Convolution Kernel Analysis (for Gridding)');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Clear Naming
%---------------------------------------------
inds = strcmp('Analysis_Name',{SCRPTipt.labelstr});
indnum = find(inds==1);
if length(indnum) > 1
    indnum = indnum(SCRPTGBL.RWSUI.scrptnum);
end
SCRPTipt(indnum).entrystr = '';
setfunc = 1;
DispScriptParam(SCRPTipt,setfunc,SCRPTGBL.RWSUI.tab,SCRPTGBL.RWSUI.panelnum);

%---------------------------------------------
% Tests
%---------------------------------------------
if not(isfield(SCRPTGBL,'Kern_File_Data'))
    if isfield(SCRPTGBL.CurrentTree.('Kern_File').Struct,'selectedfile')
    file = SCRPTGBL.CurrentTree.('Kern_File').Struct.selectedfile;
        if not(exist(file,'file'))
            err.flag = 1;
            err.msg = '(Re) Load Kern_File - path no longer valid';
            ErrDisp(err);
            return
        else
            Status('busy','Load Kern_File');
            load(file);
            saveData.path = file;
            SCRPTGBL.('Kern_File_Data') = saveData;
        end
    else
        err.flag = 1;
        err.msg = '(Re) Load Kern_File';
        ErrDisp(err);
        return
    end
end

%---------------------------------------------
% Get Input
%---------------------------------------------
ANLZ.method = SCRPTGBL.CurrentTree.Func;
ANLZ.reps = str2double(SCRPTGBL.CurrentTree.('Reps')); 
ANLZ.anlzres = str2double(SCRPTGBL.CurrentTree.('anlzRes'));
ANLZ.anlzzf = str2double(SCRPTGBL.CurrentTree.('anlzZF'));

%---------------------------------------------
% Load Implementation
%---------------------------------------------
Kern = SCRPTGBL.Kern_File_Data.KRNprms;

%---------------------------------------------
% Perform Analysis
%---------------------------------------------
func = str2func([ANLZ.method,'_Func']);
INPUT.Kern = Kern;
[ANLZ,err] = func(ANLZ,INPUT);
if err.flag
    return
end
clear INPUT;

%--------------------------------------------
% Output to TextBox
%--------------------------------------------
global FIGOBJS
FIGOBJS.(SCRPTGBL.RWSUI.tab).Info.String = ANLZ.ExpDisp;

%--------------------------------------------
% Name
%--------------------------------------------
name = ['ANLZ_',Kern.name];

%--------------------------------------------
% Return
%--------------------------------------------
name = inputdlg('Name Analysis:','Name',1,{name});
if isempty(name)
    SCRPTGBL.RWSUI.SaveGlobal = 'no';
    return
end
ANLZ.name = name{1};

SCRPTipt(indnum).entrystr = ANLZ.name;
SCRPTGBL.RWSUI.SaveVariables = ANLZ;
SCRPTGBL.RWSUI.SaveVariableNames = 'ANLZ';
SCRPTGBL.RWSUI.SaveGlobal = 'yes';
SCRPTGBL.RWSUI.SaveGlobalNames = ANLZ.name;
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = 'outloc';
SCRPTGBL.RWSUI.SaveScriptName = ANLZ.name;

Status('done','');
Status2('done','',2);
Status2('done','',3);
