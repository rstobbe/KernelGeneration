%=========================================================
% 
%=========================================================

function [default] = Inv_KB3DCart_v1c_Default2(SCRPTPATHS)

width = '5';
beta = '11';
SS = '1.25';
ZF = '128';
anlzRes = '0.125';
anlzZF = '408';
Elip = 1;

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'IF_Name';
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
default{m,1}.labelstr = 'SubSamp';
default{m,1}.entrystr = SS;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'ZeroFill';
default{m,1}.entrystr = ZF;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Elip';
default{m,1}.entrystr = Elip;

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


