%==================================================
% (v1f) 
%     Start from 'Kern_FiltCube_v1g'
%==================================================

function [SCRPTipt,KDES,err] = Kern_Cube_v1g(SCRPTipt,KDESipt)

Status2('done','SglC_KB3DCart Kernel',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Return Input
%---------------------------------------------
KDES.method = KDESipt.Func;
KDES.CubeWid = str2double(KDESipt.('CubeWid'));
KDES.Res = str2double(KDESipt.('Res'));
KDES.SubSamp = str2double(KDESipt.('SubSamp')); 

Status2('done','',2);
Status2('done','',3);

