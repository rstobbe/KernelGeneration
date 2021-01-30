%=========================================================
% 
%=========================================================

function [default] = Inv_KB3DCartS_v1c_Default2(SCRPTPATHS)

width = '5';
beta = '11';
SS = '1.25';
ZF = '128';
anlzRes = '0.125';
anlzZF = '408';

m = 1;
default{m,1}.entrytype = 'OutputName';
default{m,1}.labelstr = 'InvFilt_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'ScriptName';
default{m,1}.labelstr = 'Script_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Width';
default{m,1}.entrystr = width;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Beta';
default{m,1}.entrystr = beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'UseWithSS';
default{m,1}.entrystr = SS;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'UseWithZF';
default{m,1}.entrystr = ZF;

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
default{m,1}.scrpttype = 'InvFiltKB3DCart';
default{m,1}.labelstr = 'Generate_InvFilt';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';


