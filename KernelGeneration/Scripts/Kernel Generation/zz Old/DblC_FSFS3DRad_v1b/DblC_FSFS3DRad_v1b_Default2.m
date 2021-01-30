%=========================================================
% 
%=========================================================

function [default] = DblC_FSFS3DRad_v1b_Default2(SCRPTPATHS)

FS1_W = '12';
FS1_BW = '1.4';
FS1_beta = '9';
FS1_Res = '0.03125';

FS2_W = '12';
FS2_BW = '1.4';
FS2_beta = '9';
FS2_Res = '0.03125';

Kern_zWadd = '0.5';
SS = '2.0';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS1_Width';
default{m,1}.entrystr = FS1_W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS1_BW';
default{m,1}.entrystr = FS1_BW;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS1_Beta';
default{m,1}.entrystr = FS1_beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS1_Res';
default{m,1}.entrystr = FS1_Res;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS2_Width';
default{m,1}.entrystr = FS2_W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS2_BW';
default{m,1}.entrystr = FS2_BW;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS2_Beta';
default{m,1}.entrystr = FS2_beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FS2_Res';
default{m,1}.entrystr = FS2_Res;

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