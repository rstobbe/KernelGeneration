function ConvKernSel

W = 4;
beta = 6.5;
order = 0;
res = 0.05;

u = (-(W/2):res:(W/2));                                              % defined over |u| = W/2;     
M = beta * sqrt(1 - (2*abs(u)/W).^2);
KB = besseli(order,M);

%figure(1)
%plot(u,KB,'-*');

%KB = KB.*exp(-1i*2*pi*u*1);
KB = [zeros(1,1000-(length(KB)-1)/2) KB zeros(1,1000-(length(KB)-1)/2)];

kb = fftshift(ifftn(ifftshift(KB)));
kb = kb/max(kb(:));

t = (-2:1/100:2)*1.2;
%t = (-2:1/100:2)*2.0;
figure(2)
hold on
plot(t,real(kb(801:1201)),'b')
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([0.25 0.25],[-1 2],'k:');
plot([-0.25 -0.25],[-1 2],'k:');
plot([0.125 0.125],[-1 2],'k:');
plot([-0.125 -0.125],[-1 2],'k:');
%plot([0.6 0.6],[-1 2],'k:');
%plot([-0.6 -0.6],[-1 2],'k:');
ylabel('Conv Kernal FT Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'YLim',[-0.1 1.1]);
set(gca,'XLim',[-0.6 0.6]);
%set(gca,'XLim',[-2.4 2.4]);
%set(gcf,'units','inches');
%set(gca,'units','inches');
%set(gcf,'position',[4 4 3.2 2.8]);
%set(gca,'position',[0.5 0.5 2.5 1.8]);
%set(gcf,'paperpositionmode','auto');
%set(gca,'fontsize',10,'fontweight','bold');
%set(gca,'PlotBoxAspectRatio',[1.1 1 1]);     
%set(gca,'box','on');
