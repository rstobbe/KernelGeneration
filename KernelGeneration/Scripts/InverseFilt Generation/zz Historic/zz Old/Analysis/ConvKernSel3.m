function ConvKernSel3

W = 6;
beta = 10.5;
order = 0;
res = 0.05;

u = (-(W/2):res:(W/2));                                              % defined over |u| = W/2;     
M = beta * sqrt(1 - (2*abs(u)/W).^2);
KB = besseli(order,M);

%figure(1)
%plot(u,KB,'-*');

KB1 = [zeros(1,1000-(length(KB)-1)/2) KB zeros(1,1000-(length(KB)-1)/2)];
kb1 = fftshift(ifftn(ifftshift(KB1)));
kb0 = kb1/max(kb1(:));

kb = 0;
N = [(-8:-1) (1:8)];
for n = 1:length(N)
    KB1 = KB.*exp(-1i*2*pi*u*N(n));
    KB1 = [zeros(1,1000-(length(KB)-1)/2) KB1 zeros(1,1000-(length(KB)-1)/2)];
    kb1 = fftshift(ifftn(ifftshift(KB1)));
    kb = kb1/max(kb1(:)) + kb;
end

Err = (kb)./kb0;

t = (-2:1/100:2)*1.2;
figure(3)
hold on
%plot(t,real(Err(801:1201)),'k')
%plot(t,real(kb0(801:1201)),'k')
plot(t,real(kb(801:1201)),'k')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([0.6 0.6],[-1 2],'k:');
plot([-0.6 -0.6],[-1 2],'k:');
%plot(imag(kb))
ylabel('Conv Kernal Err','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
%set(gca,'YLim',[0.95 1.05]);
set(gca,'YLim',[-0.05 0.05]);
%set(gca,'XLim',[-2.4 2.4]);
set(gca,'XLim',[-0.7 0.7]);
set(gcf,'units','inches');
set(gca,'units','inches');
set(gcf,'position',[4 4 3.2 2.8]);
set(gca,'position',[0.5 0.5 2.5 1.8]);
set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',10,'fontweight','bold');
set(gca,'PlotBoxAspectRatio',[1.1 1 1]);     
set(gca,'box','on');
