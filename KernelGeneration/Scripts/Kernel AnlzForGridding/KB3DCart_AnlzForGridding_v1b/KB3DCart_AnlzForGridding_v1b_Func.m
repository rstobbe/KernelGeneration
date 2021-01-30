%==================================================
% (v1a)
%
%==================================================

function [ANLZ,err] = KB3DCart_AnlzForGridding_v1b_Func(ANLZ,INPUT)

Status('busy','Convolution Kernel Analysis (for Gridding)');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Get Input
%---------------------------------------------
Kern = INPUT.Kern;
clear INPUT

clr = 'r';

%---------------------------------------------
% Get Info
%---------------------------------------------
W = Kern.W;
beta = Kern.beta;
SS = Kern.DesforSS; 
Reps = ANLZ.reps;
anlzres = ANLZ.anlzres;
anlzzf = ANLZ.anlzzf;

if rem(round(1e9*(W*SS/2))/1e9,anlzres)
    err.flag = 1;
    err.msg = 'Width*TestForSS/2 not a multiple of anlzRes';
    return
end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',2);
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

figure(1000); subplot(2,2,1); hold on;
plot(t,real(kb0),clr)
plot(t,imag(kb0),'b')
plot([0.5 0.5],[-0.2 1.2],'k:');
plot([-0.5 -0.5],[-0.2 1.2],'k:');
plot([SS/2 SS/2],[-0.2 1.2],'k:');
plot([-SS/2 -SS/2],[-0.2 1.2],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);

%-------------------------------------------------
% Figure size
%-------------------------------------------------
h = gcf;
h.Units = 'Inches';
h.Position = [5 1 10 8];

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
    
    figure(1000); subplot(2,2,2); hold on;
    plot(t,real(kb1),clr);
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
fh = figure(1000); 
subplot(2,2,3); hold on;
plot(t(ind1:ind2),real(Err(ind1:ind2)),clr)
plot(t(ind1:ind2),imag(Err(ind1:ind2)),'b')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([SS/2 SS/2],[-1 2],'k:');
plot([-SS/2 -SS/2],[-1 2],'k:');
ylabel('Relative Conv Kernal Err','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);
set(gca,'YLim',[-0.05 0.05]);

subplot(2,2,4); hold on;
plot(t(ind1:ind2),real(Err(ind1:ind2)),clr)
plot(t(ind1:ind2),imag(Err(ind1:ind2)),'b')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([SS/2 SS/2],[-1 2],'k:');
plot([-SS/2 -SS/2],[-1 2],'k:');
ylabel('Relative Conv Kernal Err','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);
%set(gca,'YLim',[-0.05 0.05]);

ANLZ.Figure(1).Name = 'Kernel Characteristics';
ANLZ.Figure(1).Type = 'Graph';
ANLZ.Figure(1).hFig = fh;
ANLZ.Figure(1).hAx = gca;

%--------------------------------------
% Return
%--------------------------------------
ANLZ.ExpDisp = '';

Status2('done','',1);
Status2('done','',2);
Status2('done','',3);



