%=========================================================
% 
%=========================================================

function [default] = KB3DCart_AnlzForGridding_v1a_Default2(SCRPTPATHS)

Reps = '12';
anlzRes = '0.0125';
anlzZF = '64000';

m = 1;
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

