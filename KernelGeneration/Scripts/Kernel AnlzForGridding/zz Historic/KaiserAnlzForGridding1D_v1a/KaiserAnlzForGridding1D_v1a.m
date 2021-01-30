%==================================================
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = KaiserAnlzForGridding1D_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

KernWidth = str2double(SCRPTipt(strcmp('KernWidth',{SCRPTipt.labelstr})).entrystr);
KernBeta = str2double(SCRPTipt(strcmp('KernBeta',{SCRPTipt.labelstr})).entrystr); 
KernRes = str2double(SCRPTipt(strcmp('KernRes',{SCRPTipt.labelstr})).entrystr); 
SS = str2double(SCRPTipt(strcmp('SubSamp',{SCRPTipt.labelstr})).entrystr); 
ZF = str2double(SCRPTipt(strcmp('ZF',{SCRPTipt.labelstr})).entrystr); 
res1 = str2double(SCRPTipt(strcmp('Res_dev',{SCRPTipt.labelstr})).entrystr); 
zf1 = str2double(SCRPTipt(strcmp('ZF_dev',{SCRPTipt.labelstr})).entrystr); 

if rem(ZF,2)
    err.flag = 1;
    err.msg = 'ZF must be even';
    return
end
if rem((KernWidth/2),KernRes)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(round(1e9*(1/(KernRes*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end
if rem(round(1e9*(KernWidth*SS/2))/1e9,res1)
    err.flag = 1;
    err.msg = 'KernWidth*SS/2 not a multiple of Res_dev';
    return
end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',1);
W = KernWidth*SS;
u = abs(-(W/2):res1:(W/2)-res1);                                                  
M = KernBeta * sqrt(1 - (2*u/W).^2);
KB = besseli(0,M);

FOV = SS/(2*res1);                                                        
t = (-FOV:(SS/res1)/zf1:FOV-(SS/res1)/zf1);

KB1 = [zeros(1,zf1/2-(length(KB)/2)) KB zeros(1,zf1/2-(length(KB)/2))];
kb1 = fftshift(ifftn(ifftshift(KB1)));
kb0 = kb1/max(abs(kb1(:)));

%figure(1); hold on;
%plot(t,real(kb0),'r');
%plot(t,imag(kb0),'b');

N = [(-8:-1) (1:8)];            % 8 replications
kb = 0;
for n = 1:length(N)
    KB1 = KB.*exp(-1i*2*pi*u*N(n));
    KB1 = [zeros(1,zf1/2-(length(KB)/2)) KB1 zeros(1,zf1/2-(length(KB)/2))];
    kb1 = fftshift(ifftn(ifftshift(KB1)));
    %figure(2); hold on;
    %plot(t,real(kb1),'r');
    %plot(t,imag(kb1),'b');
    
    kb = kb1/max(abs(kb1(:))) + kb;
    %figure(3); clf; hold on
    %plot(t,real(kb),'r');
    %plot(t,imag(kb),'b');
end

Err = (kb)./kb0;

ind1 = find(t>=-SS/2,1,'first');
ind2 = find(t>=SS/2,1,'first');

figure(3)
hold on
plot(t(ind1:ind2),real(Err(ind1:ind2)),'r')
%plot(t(ind1:ind2),imag(Err(ind1:ind2)),'b')
%plot(t(ind1:ind2),real(kb0(ind1:ind2)),'k')
%plot(t(ind1:ind2),real(kb(ind1:ind2)),'k')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([1 1],[-1 2],'k:');
plot([-1 -1],[-1 2],'k:');
plot([SS/2 SS/2],[-1 2],'k:');
plot([-SS/2 -SS/2],[-1 2],'k:');
%plot(imag(kb))
ylabel('Conv Kernal Err','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
%set(gca,'YLim',[-0.05 0.05]);
%set(gca,'YLim',[-0.1 0.1]);
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);
%set(gcf,'units','inches');
%set(gca,'units','inches');
%set(gcf,'position',[4 4 3.2 2.8]);
%set(gca,'position',[0.5 0.5 2.5 1.8]);
%set(gcf,'paperpositionmode','auto');
%set(gca,'fontsize',10,'fontweight','bold');
%set(gca,'PlotBoxAspectRatio',[1.1 1 1]);     
%set(gca,'box','on');

