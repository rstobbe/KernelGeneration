

function Analysis_IF

path = 'D:\8 Programs\2 NL-PRODS\Inversion Filters\Current Files';
[file,path] = uigetfile(path);
load(strcat(path,file));

t = (1:ZF);

%figure(4);
%plot(t,10*log(V(:,ZF/2+1,ZF/2+1)));
%hold on
%title('Log Convolution Roll-off Across Full Image Space');
%drawnow;

figure(5);
if rem(ZF+1,2)
    plot(t,V(:,ZF/2+1,ZF/2+1),'r');
else
    plot(t,V(:,(ZF+1)/2,(ZF+1)/2),'r');
end
set(gca,'ylim',[0 1]);
hold on
title('Convolution Roll-off Across Full Image Space');
drawnow;


