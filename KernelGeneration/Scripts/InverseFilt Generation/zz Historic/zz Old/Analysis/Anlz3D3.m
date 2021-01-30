%=========================================================
% Generate Inverse Filter for Kaiser_Bessel Regridding
%=========================================================

function Anlz3D3

global kb

%==================
W = 6;                  
beta = 18;              
order = 0;
%==================

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
res = 0.1;
KB = KaiBesImg_Full(W,beta,res,order);
np = W/res + 1;

u = (-(W/2):res:(W/2)); 
[X,Y,Z] = meshgrid(u,u,u);

kb = 0;
%N = [(-4:-1) (1:4)];
%N = [-2 -1 0 1 2];
%M = [-2 -1 0 1 2];
%P = [-2 -1 0 1 2];
N = [-1 0 1];
M = [-1 0 1];
P = [-1 0 1];
i = 0;
for n = 1:length(N)
    for m = 1:length(M)
        for p = 1:length(P)
            if N(n)==0 && M(m)==0 && P(p)==0
                test = 0
            else
                KB1 = KB.*exp(-j*2*pi*(Y*N(n)+X*M(m)+Z*P(p)));
                zf1 = 120;
                KB1 = ifftshift(KB1);                                                 % ifftshift because odd (gives max on left)
                V = Kzerofill_isotropic(KB1,zf1);
                V = fftshift(fftn(V));                              
                kb = V/max(V(:)) + kb;
                i = i+1
            end
        end
    end
end

figure(20);
hold on
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1)*1.2;
plot(t,real(kb(:,zf1/2+1,zf1/2+1)),'b');                                   % imaginary component should be ~ zero
plot(t,real(kb(:,1,1)),'r'); 
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

figure(21)
cmap = colormap('jet');
T = real(kb(zf1/2+1+(-5:5),zf1/2+1+(-5:5),zf1/2+1));
imshow(T,[-0.1 0.2],'colormap',cmap);