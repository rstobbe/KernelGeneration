%=========================================================
% Generate Inverse Filter for Kaiser_Bessel Regridding
%=========================================================

function Anlz3D

%==================
W = 6;                  
beta = 10.5;              
order = 0;
ZF = 128;
%==================

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
res = 0.1;
KB = KaiBesImg_Full(W,beta,res,order);
np = W/res + 1;

u = (-(W/2):res:(W/2)); 

[X,Y,Z] = meshgrid(u,u,u);
KB = KB.*exp(-j*2*pi*(Y*0 + X*0 + Z*0));

%-------------------------------------------------
% Shape of Filter in Image-Space (full width)
%-------------------------------------------------  
zf1 = 120;
KB = ifftshift(KB);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(KB,zf1);
V = fftshift(fftn(V));                              
V = V/max(V(:));

figure(2);
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1)*1.2;
plot(t,real(V(:,zf1/2+1,zf1/2+1)));                                   % imaginary component should be ~ zero
%plot(t,real(V(zf1/2+1,:,zf1/2+1))); 
%plot(t,real(squeeze(V(zf1/2+1,zf1/2+1,:))));
hold on
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

