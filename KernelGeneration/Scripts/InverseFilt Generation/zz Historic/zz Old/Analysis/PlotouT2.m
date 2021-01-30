function PlotouT2

global kb

%==================
W = 6;                  
beta = 10.5;              
order = 0;
res = 0.1;
zf1 = 120;
%==================

KB = KaiBesImg_Full(W,beta,res,order);
KB1 = ifftshift(KB);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(KB1,zf1);
V = fftshift(fftn(V));                              
V = V/max(V(:));
V = kb./V;

figure(21)
colormap('jet');
T = real(V(zf1/2+1+(-6:6),zf1/2+1+(-6:6),zf1/2+1+(0)));
h = image(T);
colorbar
set(h,'CDataMapping','scaled');
caxis([-0.05 0.05]);
colorbar;  

res = 0.1;
figure(20);
hold on
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1)*1.2;
plot(t,real(V(:,zf1/2+1,zf1/2+1)),'b');                                   % imaginary component should be ~ zero
plot([0.5 0.5],[-1 2],'k:');
plot([-0.5 -0.5],[-1 2],'k:');
plot([0.6 0.6],[-1 2],'k:');
plot([-0.6 -0.6],[-1 2],'k:');
ylabel('Conv Kernal FT Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'YLim',[-0.1 1.1]);
set(gca,'XLim',[-2.4 2.4]);
set(gcf,'units','inches');
set(gca,'units','inches');
set(gcf,'position',[4 4 3.2 2.8]);
set(gca,'position',[0.5 0.5 2.5 1.8]);
set(gcf,'paperpositionmode','auto');
set(gca,'fontsize',10,'fontweight','bold');
set(gca,'PlotBoxAspectRatio',[1.1 1 1]);     
set(gca,'box','on');