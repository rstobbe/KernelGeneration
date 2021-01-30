%=========================================================
% 
%=========================================================

function [default] = SglC_KB3DCart_v1b_Default2(SCRPTPATHS)

width = '5';
beta = '11';
SS = '1.25';
res = '0.01';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Width';
default{m,1}.entrystr = width;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Beta';
default{m,1}.entrystr = beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'DesForSS';
default{m,1}.entrystr = SS;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Res';
default{m,1}.entrystr = res;

m = m+1;
default{m,1}.entrytype = 'RunScrptFunc';
default{m,1}.scrpttype = 'KB3DCartGen';
default{m,1}.labelstr = 'Generate_Kernel';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';

