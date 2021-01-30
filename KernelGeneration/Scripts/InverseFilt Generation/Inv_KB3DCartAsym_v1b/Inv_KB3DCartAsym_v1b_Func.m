%==================================================
%
%==================================================

function [OUTPUT,err] = Inv_KB3DCart_v1c_Func(INPUT)

Status('busy','Create Inverse Filter');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';
OUTPUT = struct();

%---------------------------------------------
% Get Input
%---------------------------------------------
IF = INPUT.IF;
clear INPUT;

%---------------------------------------------
% Common Variables
%---------------------------------------------
W = IF.W;
beta = IF.beta;
SS = IF.SS;
ZF = IF.ZF;
Elip = IF.Elip;
anlzres = IF.anlzres;
anlzzf = IF.anlzzf;

%---------------------------------------------
% Tests
%---------------------------------------------
if rem(ZF,2)
    err.flag = 1;
    err.msg = 'ZF must be even';
    return
end
if rem(round(1e9*(W*SS/2))/1e9,anlzres)
    err.flag = 1;
    err.msg = 'Width*UseWithSS/2 not a multiple of anlzRes';
    return
end
if rem((anlzzf*anlzres)+1,2)
    err.flag = 1;
    err.msg = 'anlzZF*anlzRes must be odd';
    return
end
if (anlzzf*anlzres) > ZF
    err.flag = 1;
    err.msg = 'anlzZF*anlzRes greater than UseWithZF';
    return
end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',1);
ssW = W*SS;
u = (-(ssW/2):anlzres:(ssW/2)-anlzres);                                                  
M = beta * sqrt(1 - (2*abs(u)/ssW).^2);
KB_1D = besseli(0,M);
KB_1D = KB_1D/max(KB_1D);

N = length(KB_1D);
KB = zeros(N,N,N);
for a = 1:N
    for b = 1:N
        for c = 1:N
            KB(a,b,c) = KB_1D(a)*KB_1D(b)*KB_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2); 
end
Status2('done','',2);

figure(100); subplot(2,2,1);
cen = length(KB)/2+1;
plot(squeeze(KB(cen,cen,:)));                                   
ylabel('kSpace','fontsize',10,'fontweight','bold');
xlabel('anlzRes Sampling','fontsize',10,'fontweight','bold');
title('Quantization Test');
f = gcf;
f.Units = 'inches';
f.Position = [4 1 12 9];

%-------------------------------------------------
% Shift and Zero-Fill
%-------------------------------------------------  
KB = ifftshift(KB); 
V = zerofill_isotropic_even_doubles(KB,anlzzf);

%figure(2);
%plot(squeeze(V(1,1,:)));                                   
%drawnow;

%-------------------------------------------------
% FT
%-------------------------------------------------  
Status2('busy','FT',1);
V = fftshift(ifftn(V));                              
V = V/max(max(max(V)));
FOV = SS/(2*anlzres);                                                       
t = (-FOV:(SS/anlzres)/anlzzf:FOV-(SS/anlzres)/anlzzf);

figure(100); subplot(2,2,2); hold on;
plot(t,imag(V(:,anlzzf/2+1,anlzzf/2+1)),'b');                                 % imaginary component should be ~ zero
plot(t,real(V(:,anlzzf/2+1,anlzzf/2+1)),'r');                                 % real component - 'minimal fold-over at edges'
plot(t,squeeze(real(V(anlzzf/2+1,:,anlzzf/2+1))),'r'); 
plot(t,squeeze(real(V(anlzzf/2+1,anlzzf/2+1,:))),'r'); 
ylim([-0.1 1]);
ylabel('ImageSpace Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
title('Fold Test');

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = (anlzzf/2)*(1-anlzres)+1.5;
c = (anlzzf/2)*(1+anlzres)+0.5;
V2 = V(b:c,b:c,b:c);
V2 = V2/max(max(max(V2)));
np2 = length(V2);

figure(100); subplot(2,2,3); hold on;
t2 = t(b:c);
plot(t2,real(V2(:,(np2+1)/2,(np2+1)/2)),'r-*');
ylabel('ImageSpace Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
title('Quantization Test');

%-------------------------------------------------
% Test Filter Inversion 
%-------------------------------------------------
Image = 1./V2;
bot = length(Image)*(SS-1)/(2*SS)+0.5;
top = length(Image)*(SS+1)/(2*SS)-0.5;
Image = Image(bot:top,bot:top,bot:top);

top = 100;
MSTRCT.dispwid = [0 top];
INPUT.Image = Image;
INPUT.MSTRCT = MSTRCT;
PlotMontage_v1c(INPUT);

%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
Status2('busy','Interpolate to ZF',1);
V2 = fftshift(fftn(ifftshift(V2)));

%figure(5);
%np3 = length(V2);
%plot(1:np3,real(V2(:,(np3+1)/2,(np3+1)/2)),'r-');
%hold on
%title('Kernel inter');
%drawnow;

%V = zerofill_isotropic_odd_doubles(ifftshift(V2),ZF);
for n = 2:2:20
    if not(rem(n/SS,1))
        step = n;
        break
    end
end
OutZF = step*round(ZF*Elip/step);
Elip = OutZF/ZF;
zfdims = [ZF,ZF,OutZF];
V = zeros(zfdims);
sz = size(V2);
bot = floor((zfdims - sz(1:3))/2)+1;
top = bot + sz(1:3) - 1;
V(bot(1):top(1),bot(2):top(2),bot(3):top(3)) = V2;

%figure(6); hold on
%plot(1:ZF,real(V(:,1,1)),'r-*');
%hold on
%title('Zerofill Kernel');
%drawnow;

V = fftshift(ifftn(ifftshift(V)));
V = V/max(max(max(V)));
V = abs(V);

figure(100); subplot(2,2,3); hold on;
t = (-0.5*SS:SS/ZF:0.5*SS-SS/ZF);
plot(t,V(:,ZF/2+1,OutZF/2+1),'c');
plot([-0.5 -0.5],[-0.2 1],'k:'); 
plot([0.5 0.5],[-0.2 1],'k:'); 
plot([-0.5*SS 0.5*SS],[0 0],'k:');
axis([-0.5*SS 0.5*SS -0.1 1]);

%-------------------------------------------------
% Output
%-------------------------------------------------
IFprms = IF;
IFprms.W = W;
IFprms.beta = beta;
IFprms.SS = SS;
IFprms.ZF = ZF;
IFprms.Elip = Elip;
IFprms.anlzres = anlzres;
IFprms.anlzzf = anlzzf;
IFprms.V = V;
name = ['IF_KBCw',num2str(W,2),'b',num2str(beta,3),'ss',num2str(SS,3),'zf',num2str(ZF,3),'e',num2str(Elip,2)];
IFprms.name = regexprep(name,'\.', 'p');

OUTPUT.IFprms = IFprms;

