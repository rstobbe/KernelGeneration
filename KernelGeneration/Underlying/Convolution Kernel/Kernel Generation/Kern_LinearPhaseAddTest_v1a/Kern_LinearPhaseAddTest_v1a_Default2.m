%=========================================================
% 
%=========================================================

function [default] = Kern_PhaseAddTest_v1a_Default2(SCRPTPATHS)

W = '12';
BW = '1.4';
beta = '9';
Res = '0.0125';
SS = '2.0';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'Width';
default{m,1}.entrystr = W;

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
default{m,1}.labelstr = 'Res';
default{m,1}.entrystr = Res;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'DesForSS';
default{m,1}.entrystr = SS;

