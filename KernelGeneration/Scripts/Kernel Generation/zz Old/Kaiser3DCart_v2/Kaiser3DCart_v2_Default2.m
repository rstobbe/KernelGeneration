%=========================================================
% 
%=========================================================

function [default] = Kaiser3DCart_v2_Default2(SCRPTPATHS)

kernwid = '5';
kernres = '0.01';
beta = '11';
zW_add = 1.5;
SS = 1.25;

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