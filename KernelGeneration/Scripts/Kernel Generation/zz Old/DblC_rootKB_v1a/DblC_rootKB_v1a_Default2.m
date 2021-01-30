%=========================================================
% 
%=========================================================

function [default] = DblC_rootKB_v1a_Default2(SCRPTPATHS)

KB_W = '6';
KB_beta = '12';
KB_Res = '0.0075';

Kern_zWadd = '0.6';
SS = '1.25';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KB_Width';
default{m,1}.entrystr = KB_W;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KB_Beta';
default{m,1}.entrystr = KB_beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KB_Res';
default{m,1}.entrystr = KB_Res;

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