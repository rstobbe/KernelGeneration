%=========================================================
% 
%=========================================================

function [KRNprms,err] = DblC_FSsq3DRad_v1f_Func(INPUT,KRNprms)

Status('busy','Create Kernel');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Get Input
%---------------------------------------------
W = KRNprms.W;
BW = KRNprms.BW;
res = KRNprms.res;
beta = KRNprms.beta;
SS = KRNprms.DesforSS;
clear INPUT;

%---------------------------------------------
% Tests
%---------------------------------------------
if rem((W/2),res)
    err.flag = 1;
    err.msg = 'Width/2 not a multiple of Res';
    return
end
if rem(round(1e9*(1/(res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(Res*DesForSS) not an integer';
    return
end

%---------------------------------------------
% Find all SS values that are compatible with res
%---------------------------------------------
m = 1;
for n = 1:0.0005:6
    test = round(1e9/(res*n))/1e9;
    if not(rem(test,1)) && not(rem(W*n,1))
        posSS(m) = n;
        m = m+1;
    end    
end

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
for n = 1:length(posSS)
    zWtest = 2*(ceil((W*posSS(n)-2)/2)+1)/posSS(n);
    zW = W + 0.1;
    while true
        if zW > zWtest
            break
        end
        zW = zW + 0.1;
    end
    zWarr(n) = round(zW*1e9)/1e9;
end
zW = max(zWarr);

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
if (ceil((zW/2)*(1/res))+1) > 1100
    err.flag = 1;
    err.msg = ['kernel size: ',num2str(ceil((zW/2)*(1/res))+1)];
    return
end

%--------------------------------------------
% Generate FS Kernel Shape
%--------------------------------------------
res0 = 0.0001;
u0 = (0:res0:(W/2));   
FTSphereProf1 = 3*(sin(u0*pi*BW) - pi*BW*u0.*cos(u0*pi*BW))./((u0*pi*BW).^3);
FTSphereProf1(1) = 1;                                              
M = beta * sqrt(1 - (2*u0/W).^2);
KBProf1 = besseli(0,M);
KBProf1 = KBProf1/max(KBProf1);
KernProf1 = FTSphereProf1.*KBProf1;

figure(200); hold on;
%plot(u0,FTSphereProf1,'m');
%plot(u0,KBProf1,'g');
plot(u0,KernProf1,'b');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

%--------------------------------------------
% CSV Single Kernel
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value of Single Kernel',1);
CSVres = 0.05;
%CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5; 
uT = (-(W/2):CSVres:(W/2)); 
Rad = ones(length(uT),length(uT),length(uT));
tKern1 = zeros(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                tKern1(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
FwdKern_CSV = sum(tKern1(:))*CSVres^3;
RvsKern_CSV = FwdKern_CSV;

interptest = 0;
if interptest == 1
    figure(200); plot(uT,tKern1(:,(N+1)/2,(N+1)/2),'k*');
    %figure(200); plot(Rad(:)*(W/2),tKern1(:),'k*');
    %figure(200); plot(Rad(:)*(W/2),tKern2(:),'k*');
end

%--------------------------------------------
% Determine Double Kern k-Space
%--------------------------------------------
Status2('busy','Determine Double Kern k-Space Shape',1);
ZF = 320;
C = ZF/2 + 1;
TF1 = fftshift(fftn(ifftshift(tKern1)));                                     
TF = TF1.^2;
ZFTF = zerofill_isotropic_odd_doubles(ifftshift(TF),ZF);
DblKern = fftshift(ifftn(ZFTF));
clear ZFTF
DblKernProf = DblKern(:,C,C)/abs(DblKern(C,C,C));
clear DblKern
iW = W + CSVres;
tu0 = (-(iW/2):iW/ZF:(iW/2)-iW/ZF);
iDblKernProf = interp1(tu0,DblKernProf,u0);

figure(200); hold on;
plot(u0,iDblKernProf,'r');
DblKernProf = real(iDblKernProf);

%--------------------------------------------
% Determine Double Kern Shape
%--------------------------------------------
Status2('busy','Determine Double Kern Image Shape',1);
ZFtKern1 = zerofill_isotropic_odd_doubles(ifftshift(tKern1),ZF);
TF1 = fftshift(fftn(ZFtKern1));
clear ZFtKern1
TF = TF1.^2;
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF1
iCSVres = 1/((1/CSVres)+(1/W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));

figure(300); hold on;
plot(f,real(TF1Prof),'b');
plot(f,imag(TF1Prof),'b:');
plot(f,real(TFProf),'r');
plot(f,imag(TFProf),'r:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

imageshow = 0;
if imageshow == 1
    tTF = TF(ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6));
    rows = floor(sqrt(length(tTF)));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1;  IMSTRCT.stop = length(tTF); 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 max(abs(tTF(:)))*1.05]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; 
    IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 400; 
    AxialMontage_v2a(tTF,IMSTRCT);
end

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value of Double Kernel',1);
Rad = ones(length(uT),length(uT),length(uT));
tKern = zeros(length(uT),length(uT),length(uT));
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                tKern(a,b,c) = lin_interp4(DblKernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
DblKern_CSV = sum(tKern(:))*CSVres^3;

%--------------------------------------------
% Test
%--------------------------------------------
Status2('busy','Test',1);
ZFtKern = zerofill_isotropic_odd_doubles(ifftshift(tKern),ZF);
TF = fftshift(fftn(ZFtKern));
clear ZFtKern
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF
figure(300); hold on; title('Image Profile');
plot(f,real(TFProf),'k:');
ind1 = find(f<-2,1,'last');
ind2 = find(f>2,1,'first');
xlim([f(ind1) f(ind2)]);
clear tKern
clear tKern1

Status2('done','',2);
button = questdlg('Continue?','Convolution Kernel Generation','Yes','No','Yes');
if strcmp(button,'No')
    return
end

%==============================================================================

%--------------------------------------------
% Generate Double Convolution Kernel Approximation
%--------------------------------------------
% Status2('busy','Generate Double Convolution Kernel Approximation',1);
% u = (0:res:(W/2)); 
% DblKern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
% Rad = ones(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
% N = length(u);
% for a = 1:N
%     for b = 1:N
%         for c = 1:N
%             rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
%             if rad <= 1
%                 Rad(a,b,c) = rad;
%                 DblKern(a,b,c) = lin_interp4(DblKernProf,rad,L);
%             end
%         end
%     end
%     Status2('busy',num2str(N-a),2);    
% end
% Status2('done','',2);
% interptest = 0;
% if interptest == 1
%     figure(200); plot(Rad(:)*(W/2),DblKern(:),'k*');
% end

%--------------------------------------------
% Generate FS1 Kernel
%--------------------------------------------
Status2('busy','Generate FS1 Kernel',1);
u = (0:res:(W/2)); 
FS1Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
Rad = ones(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                FS1Kern(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),FS1Kern(:),'k*');
end

%--------------------------------------------
% Output
%--------------------------------------------
KRNprms.type = 'DblC_FSsq3DRad_v1f';
KRNprms.prec = 'singles';
KRNprms.posSS = posSS;

KRNprms.DblKern.convscaleval = DblKern_CSV;
KRNprms.DblKern.OKforSS = SS;
KRNprms.DblKern.W = W;
KRNprms.DblKern.iKern = 1/res;
KRNprms.DblKern.res = res;
KRNprms.DblKern.zW = zW;
%KRNprms.DblKern.Kern = DblKern;

KRNprms.FwdKern.convscaleval = FwdKern_CSV;
KRNprms.FwdKern.W = W;
KRNprms.FwdKern.BW = BW;
KRNprms.FwdKern.beta = beta;
KRNprms.FwdKern.iKern = 1/res;
KRNprms.FwdKern.res = res;
KRNprms.FwdKern.zW = zW;
KRNprms.FwdKern.Kern = FS1Kern;

KRNprms.RvsKern.convscaleval = RvsKern_CSV;
KRNprms.RvsKern.W = W;
KRNprms.RvsKern.BW = BW;
KRNprms.RvsKern.beta = beta;
KRNprms.RvsKern.iKern = 1/res;
KRNprms.RvsKern.res = res;
KRNprms.RvsKern.zW = zW;
%KRNprms.RvsKern.Kern = FS1Kern;   
  
name = ['Kern_FSsqRw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3),'ss',num2str(SS,3),'r',num2str(res,3)];
KRNprms.name = regexprep(name,'\.', 'p');

Status('done','');
Status2('done','',2);
Status2('done','',3);
