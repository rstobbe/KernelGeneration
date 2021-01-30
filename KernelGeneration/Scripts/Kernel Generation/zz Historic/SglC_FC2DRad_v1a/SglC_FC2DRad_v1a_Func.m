%=========================================================
% 
%=========================================================

function [KRNprms,err] = SglC_FC2DRad_v1a_Func(INPUT,KRNprms)

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
for n = 1:0.0005:10
    test = round(1e9/(res*n))/1e9;
    if not(rem(test,1))
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
kernelsize = ceil((zW/2)*(1/res))+1;
button = questdlg(['Kernel Dimensions: ',num2str(kernelsize),'.  Continue?']);
if strcmp(button,'No')
    err.flag = 4;
    err.msg = '';
    return
end

%--------------------------------------------
% Generate FC Kernel Shape
%--------------------------------------------
res0 = 0.0001;
u0 = (0:res0:(W/2));   
FTCircleProf = 2*(besselj(1,u0*pi*BW))./(u0*pi*BW);
FTCircleProf(1) = 1;                                              
M = beta * sqrt(1 - (2*u0/W).^2);
KBProf = besseli(0,M);
KBProf = KBProf/max(KBProf);
KernProf = FTCircleProf.*KBProf;

figure(200); hold on;
%plot(u0,FTCircleProf,'m');
%plot(u0,KBProf,'g');
plot(u0,KernProf,'b');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold'); 

%--------------------------------------------
% CSV
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value',1);
CSVres = 0.1;
uT = (-(W/2):CSVres:(W/2)); 
Rad = ones(length(uT),length(uT));
tKern = zeros(length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        rad = ((uT(a).^2 + uT(b).^2).^(0.5))/(W/2);
        if rad <= 1
            Rad(a,b) = rad;
            tKern(a,b) = lin_interp4(KernProf,rad,L);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
CSV = sum(tKern(:))*CSVres^2;

interptest = 0;
if interptest == 1
    figure(200); plot(uT,tKern(:,(N+1)/2),'k*');
end

%--------------------------------------------
% Determine Kern Shape
%--------------------------------------------
Status2('busy','Determine Kern Image Shape',1);
ZF = 320;
C = ZF/2 + 1;
ZFtKern = zerofill_isotropic2D_odd_doubles(ifftshift(tKern),ZF);
TF = fftshift(fftn(ZFtKern));
clear ZFtKern
TFProf = squeeze(TF(C,:))/abs(TF(C,C));
clear TF
iCSVres = 1/((1/CSVres)+(1/W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));

figure(300); hold on;
plot(f,real(TFProf),'b');
plot(f,imag(TFProf),'b:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

imageshow = 0;
if imageshow == 1
    tTF = TF;
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1;  IMSTRCT.stop = 1; 
    IMSTRCT.rows = 1; IMSTRCT.lvl = [0 max(abs(tTF(:)))*1.05]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; 
    IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 400; 
    AxialMontage_v2a(tTF,IMSTRCT);
end

Status2('done','',2);
button = questdlg('Continue?','Convolution Kernel Generation','Yes','No','Yes');
if strcmp(button,'No')
    return
end

%==============================================================================

%--------------------------------------------
% Generate Kernel 
%--------------------------------------------
Status2('busy','Generate Kernel',1);
u = (0:res:(W/2)); 
Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
Rad = ones(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        rad = ((u(a).^2 + u(b).^2).^(0.5))/(W/2);
        if rad <= 1
            Rad(a,b) = rad;
            Kern(a,b) = lin_interp4(KernProf,rad,L);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 1;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),Kern(:),'k*');
end

%----------------------------------------------------
% Panel Output
%----------------------------------------------------
Panel(1,:) = {'ConvScalVal',CSV,'Output'};
PanelOutput = cell2struct(Panel,{'label','value','type'},2);
KRNprms.PanelOutput = PanelOutput;

%--------------------------------------------
% Return
%--------------------------------------------
KRNprms.convscaleval = CSV;
KRNprms.DesforSS = SS;
KRNprms.posSS = posSS;
KRNprms.iKern = 1/res;
KRNprms.zW = zW;
KRNprms.Kern = Kern;
name = ['Kern_FCRw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3)];
KRNprms.name = regexprep(name,'\.', 'p');

Status('done','');
Status2('done','',2);
Status2('done','',3);
