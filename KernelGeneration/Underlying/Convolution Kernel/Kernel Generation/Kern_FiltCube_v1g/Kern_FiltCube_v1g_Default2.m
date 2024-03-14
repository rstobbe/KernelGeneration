%=========================================================
% 
%=========================================================

function [default] = Kern_FiltCube_v1g_Default2(SCRPTPATHS)

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'CubeWid';
default{m,1}.entrystr = 1;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FiltWid';
default{m,1}.entrystr = 2;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'FiltBeta';
default{m,1}.entrystr = 5.5;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'SubSamp';
default{m,1}.entrystr = 1.6;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Res';
default{m,1}.entrystr = 0.00625;


