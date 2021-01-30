%=========================================================
% 
%=========================================================

function [default] = KB3DCart_AnlzForGridding_v1b_Default2(SCRPTPATHS)

Reps = '12';
anlzRes = '0.0125';
anlzZF = '64000';

m = 1;
default{m,1}.entrytype = 'OutputName';
default{m,1}.labelstr = 'Analysis_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'ScriptName';
default{m,1}.labelstr = 'Script_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'RunExtFunc';
default{m,1}.labelstr = 'Kern_File';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Select';
default{m,1}.runfunc1 = 'LoadConvKernCur';
default{m,1}.(default{m,1}.runfunc1).curloc = SCRPTPATHS.outloc;
default{m,1}.runfunc2 = 'LoadConvKernDef';
default{m,1}.(default{m,1}.runfunc2).defloc = SCRPTPATHS.outloc;
default{m,1}.searchpath = SCRPTPATHS.scrptshareloc;
default{m,1}.path = SCRPTPATHS.scrptshareloc;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Reps';
default{m,1}.entrystr = Reps;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'anlzRes';
default{m,1}.entrystr = anlzRes;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'anlzZF';
default{m,1}.entrystr = anlzZF;

m = m+1;
default{m,1}.entrytype = 'RunScrptFunc';
default{m,1}.scrpttype = 'ConvKernGridTest';
default{m,1}.labelstr = 'Test_Kernel';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';

