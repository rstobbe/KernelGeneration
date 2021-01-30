%=========================================================
% 
%=========================================================

function [default] = DblC_KBsq3DRad_v1a_Default2(SCRPTPATHS)

W = '5';
beta = '5';
Res = '0.01';

zWadd = '0.6';
SS = '2.5';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Width';
default{m,1}.entrystr = W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Beta';
default{m,1}.entrystr = beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Res';
default{m,1}.entrystr = Res;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'zWadd';
default{m,1}.entrystr = zWadd;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'TestForSS';
default{m,1}.entrystr = SS;

m = m+1;
default{m,1}.entrytype = 'Choose';
default{m,1}.labelstr = 'Visuals';
default{m,1}.entrystr = 'On';
default{m,1}.options = {'Off','On'};

m = m+1;
default{m,1}.entrytype = 'RunScrptFunc';
default{m,1}.labelstr = 'Generate Kernel';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';