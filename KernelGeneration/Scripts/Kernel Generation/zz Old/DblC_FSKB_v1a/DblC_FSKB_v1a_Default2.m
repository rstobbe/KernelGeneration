%=========================================================
% 
%=========================================================

function [default] = DblC_FSKB_v1a_Default2(SCRPTPATHS)

FS_W = '10';
FS_BW = '1.1';
FS_beta = '6';

KB_W = '2';
KB_beta = '9';

Kern_Res = '0.0125';
Kern_zWadd = '0.6';
SS = '3.2';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS_Width';
default{m,1}.entrystr = FS_W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS_BW';
default{m,1}.entrystr = FS_BW;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS_Beta';
default{m,1}.entrystr = FS_beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KB_Width';
default{m,1}.entrystr = KB_W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KB_Beta';
default{m,1}.entrystr = KB_beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Kern_Res';
default{m,1}.entrystr = Kern_Res;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Kern_zWadd';
default{m,1}.entrystr = Kern_zWadd;

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