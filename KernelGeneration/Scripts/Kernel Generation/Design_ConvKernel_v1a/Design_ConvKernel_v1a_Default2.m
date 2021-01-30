%=========================================================
% 
%=========================================================

function [default] = Design_ConvKernel_v1a_Default2(SCRPTPATHS)

if strcmp(filesep,'\')
    KERNPath = [SCRPTPATHS.pioneerloc,'Other\Underlying\Convolution Kernel\Kernel Generation\'];
elseif strcmp(filesep,'/')
end
KERNfunc = 'DblC_FSsq3DRad_v2f';

m = 1;
default{m,1}.entrytype = 'OutputName';
default{m,1}.labelstr = 'Kernel_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'ScriptName';
default{m,1}.labelstr = 'Script_Name';
default{m,1}.entrystr = '';

m = m+1;
default{m,1}.entrytype = 'ScrptFunc';
default{m,1}.labelstr = 'Kernfunc';
default{m,1}.entrystr = KERNfunc;
default{m,1}.searchpath = KERNPath;
default{m,1}.path = [KERNPath,KERNfunc];

m = m+1;
default{m,1}.entrytype = 'RunScrptFunc';
default{m,1}.scrpttype = 'ConvKernel';
default{m,1}.labelstr = 'Design';
default{m,1}.entrystr = '';
default{m,1}.buttonname = 'Run';

