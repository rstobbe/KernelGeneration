%==================================================
% (v1a)
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = KB3DCart_AnlzForGridding_v1a(SCRPTipt,SCRPTGBL)

err.flag = 0;
err.msg = '';

AllTrees = SCRPTGBL.AllTrees;
if not(isfield(AllTrees,'KB3DCartGen'))
    err.flag = 1;
    err.msg = '''KB3DCartGen'' script must also be loaded';
    ErrDisp(err);
    return
end

W = str2double(SCRPTGBL.AllTrees.KB3DCartGen.Width);
beta = str2double(SCRPTGBL.AllTrees.KB3DCartGen.Beta);
SS = str2double(SCRPTGBL.AllTrees.KB3DCartGen.TestForSS); 
Reps = str2double(SCRPTGBL.CurrentScript.Reps);
anlzres = str2double(SCRPTGBL.CurrentScript.anlzRes);
anlzzf = str2double(SCRPTGBL.CurrentScript.anlzZF);

if rem(round(1e9*(W*SS/2))/1e9,anlzres)
    err.flag = 1;
    err.msg = 'Width*TestForSS/2 not a multiple of anlzRes';
    return
end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',1);
ssW = W*SS;
u = (-(ssW/2):anlzres:(ssW/2)-anlzres);                                                  
M = beta * sqrt(1 - (2*abs(u)/ssW).^2);
KB = besseli(0,M);

%-------------------------------------------------
% Profile in Image Space
%-------------------------------------------------
KB1 = [zeros(1,anlzzf/2-(length(KB)/2)) KB zeros(1,anlzzf/2-(length(KB)/2))];
kb1 = fftshift(ifftn(ifftshift(KB1)));
kb0 = kb1/max(abs(kb1(:)));

%-------------------------------------------------
% Plot Profile in Image Space
%-------------------------------------------------
FOV = SS/(2*anlzres);                                                        
t = (-FOV:(SS/anlzres)/anlzzf:FOV-(SS/anlzres)/anlzzf);
ind1 = find(t>=-SS/2,1,'first');
ind2 = find(t>=SS/2,1,'first');
figure(1201); hold on;
plot(t,real(kb0),'r')
plot(t,imag(kb0),'b')
plot([0.5 0.5],[-0.2 1.2],'k:');
plot([-0.5 -0.5],[-0.2 1.2],'k:');
plot([SS/2 SS/2],[-0.2 1.2],'k:');
plot([-SS/2 -SS/2],[-0.2 1.2],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);

%-------------------------------------------------
% Generate and Sum Replications
%-------------------------------------------------
N = [(-Reps:-1) (1:Reps)];            
kb = 0;
for n = 1:length(N)
    KB0 = KB.*exp(-1i*2*pi*u*N(n));
    KB1 = [zeros(1,anlzzf/2-(length(KB)/2)) KB0 zeros(1,anlzzf/2-(length(KB)/2))];
    kb1 = fftshift(ifftn(ifftshift(KB1)));
    kb1 = kb1/max(abs(kb1(:)));
    
    figure(1202); hold on;
    plot(t,real(kb1),'r');
    plot(t,imag(kb1),'b');
    xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
    ylabel('Kernel Replications','fontsize',10,'fontweight','bold');
    
    kb = kb1 + kb;
end

%-------------------------------------------------
% Mult by inverse of kernel
%-------------------------------------------------
Err = (kb)./kb0;

%-------------------------------------------------
% Plot Error
%-------------------------------------------------
figure(1203)
hold on
plot(t(ind1:ind2),real(Err(ind1:ind2)),'r')
plot(t(ind1:ind2),imag(Err(ind1:ind2)),'b')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([SS/2 SS/2],[-1 2],'k:');
plot([-SS/2 -SS/2],[-1 2],'k:');
ylabel('Relative Conv Kernal Err','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);
set(gca,'YLim',[-0.05 0.05]);

%figure(1204); clf; hold on
%plot(t(ind1:ind2),real(kb(ind1:ind2)),'r')
%plot(t(ind1:ind2),imag(kb(ind1:ind2)),'b')
