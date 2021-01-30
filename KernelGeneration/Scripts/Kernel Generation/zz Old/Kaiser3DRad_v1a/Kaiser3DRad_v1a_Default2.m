%=========================================================
% 
%=========================================================

function [default] = Kaiser3DRad_v1a_Default2(SCRPTPATHS)

kernwid = '10';
kernres = '0.0125';
beta = '6';
zW_add = 0.6;
SS = 2.5;

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernWidth';
default{m,1}.entrystr = kernwid;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernBeta';
default{m,1}.entrystr = beta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernRes';
default{m,1}.entrystr = kernres;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'zW_add';
default{m,1}.entrystr = zW_add;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'SubSamp';
default{m,1}.entrystr = SS;

m = m+1;
default{m,1}.entrytype = 'Choose';
default{m,1}.labelstr = 'Visuals';
default{m,1}.entrystr = 'On';
default{m,1}.options = {'Off','On'};