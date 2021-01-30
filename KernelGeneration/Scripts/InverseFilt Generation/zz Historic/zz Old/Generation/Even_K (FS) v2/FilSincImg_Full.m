%=========================================================
% 
%=========================================================

function [FS] = FilSincImg_Full(W,BW,beta,res)

if rem((W/2),res)
    error;
    return
end
if rem(1,res)
    error;
    return
end

u = abs(-(W/2):res:(W/2));                                              % defined over |u| = W/2;     
Sinc_1D = sinc(u*BW);

M = beta * sqrt(1 - (2*u/W).^2);
KB_1D = besseli(0,M);

FS_1D = Sinc_1D.*KB_1D;
FS_1D = FS_1D/max(FS_1D(:));

N = length(FS_1D);
FS = zeros(N,N,N);
for a = 1:N
    for b = 1:N
        for c = 1:N
            FS(a,b,c) = FS_1D(a)*FS_1D(b)*FS_1D(c);
        end
    end
end


