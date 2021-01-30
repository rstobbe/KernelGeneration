%///////////////////////////////////////////////////////////////////%
% 						Function: 	KaiBes Full         			%
%-------------------------------------------------------------------%
% - Create Kaiser-Bessel Filter 									%
%///////////////////////////////////////////////////////////////////%

%W = width of Kaiser-Bessel filter in grid points                                         										
%beta = window shape

function [KB] = KaiBesImg_Full(W,beta,res,order)

u = abs([-(W/2):res:(W/2)]);                                              % defined over |u| = W/2;     
M = beta * sqrt(1 - (2*u/W).^2);
KB_1D = besseli(order,M);

N = length(KB_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            KB(a,b,c) = KB_1D(a)*KB_1D(b)*KB_1D(c);
        end
    end
end

KB = KB / max(max(max(KB)));

