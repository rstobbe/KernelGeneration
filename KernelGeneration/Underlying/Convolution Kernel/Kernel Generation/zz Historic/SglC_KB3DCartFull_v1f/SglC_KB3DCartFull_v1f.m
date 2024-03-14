%==================================================
% (v1f) 
%       - Remove PossibleSS
%==================================================

function [SCRPTipt,KDES,err] = SglC_KB3DCartFull_v1f(SCRPTipt,KDESipt)

Status2('done','SglC_KB3DCart Kernel',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Return Input
%---------------------------------------------
KDES.method = KDESipt.Func;
KDES.W = str2double(KDESipt.('Width'));
KDES.res = str2double(KDESipt.('Res'));
KDES.beta = str2double(KDESipt.('Beta'));
KDES.DesforSS = str2double(KDESipt.('DesForSS')); 
KDES.prec = 'singles';

Status2('done','',2);
Status2('done','',3);

