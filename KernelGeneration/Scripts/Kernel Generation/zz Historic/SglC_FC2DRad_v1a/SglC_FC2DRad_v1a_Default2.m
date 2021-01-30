%=========================================================
% 
%=========================================================

function [default] = SglC_FS3DRad_v1e_Default2(SCRPTPATHS)

width = '10';
BW = '1';
beta = '6';
SS = 1;
res = '0.0125';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Kern_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Width';
default{m,1}.entrystr = width;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'BW';
default{m,1}.entrystr = BW;

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
default{m,1}.scrpttype = 'FS3DRadGen';
default{m,1}.labelstr = 'Generate_Kernel';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';
