%=========================================================
% 
%=========================================================

function [default] = CartKernAnlzForGridding_v1a_Default2(SCRPTPATHS)

ZF = '128';
Res_dev = '0.025';
ZF_dev = '8000';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'ZF';
default{m,1}.entrystr = ZF;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Res_dev';
default{m,1}.entrystr = Res_dev;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'ZF_dev';
default{m,1}.entrystr = ZF_dev;

m = m+1;
default{m,1}.entrytype = 'RunScrptFunc';
default{m,1}.scrpttype = 'ConvKernGridTest';
default{m,1}.labelstr = 'Test_Kernel';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';

