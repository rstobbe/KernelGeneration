%==================================================
% (v2f) -> (as v1f just lower down)
%   - Don't save double Conv
%   - allow bigger
%==================================================

function [SCRPTipt,KDES,err] = DblC_FSsq3DRad_v2f(SCRPTipt,KDESipt)

Status2('done','Get Build Kernel Design Input',2);
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
KDES.BW = str2double(KDESipt.('BW'));
KDES.DesforSS = str2double(KDESipt.('DesForSS')); 
KDES.prec = 'singles';

Status2('done','',2);
Status2('done','',3);

