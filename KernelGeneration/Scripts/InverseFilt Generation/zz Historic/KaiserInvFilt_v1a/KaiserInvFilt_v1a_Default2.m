%=========================================================
% 
%=========================================================

function [default] = KaiserInvFilt_v1a_Default2(SCRPTPATHS)

kernwid = '5';
kernbeta = '11';
kernres = '0.01';
SS = '1.25';
ZF = '128';
Res_dev = '0.125';
ZF_dev = '408';

m = 1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernWidth';
default{m,1}.entrystr = kernwid;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernBeta';
default{m,1}.entrystr = kernbeta;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'KernRes';
default{m,1}.entrystr = kernres;

m = m+1;
default{m,1}.entrytype = 'Input';
default{m,1}.labelstr = 'SubSamp';
default{m,1}.entrystr = SS;

m = m+1;
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
default{m,1}.entrytype = 'Choose';
default{m,1}.labelstr = 'Visuals';
default{m,1}.entrystr = 'On';
default{m,1}.options = {'Off','On'};