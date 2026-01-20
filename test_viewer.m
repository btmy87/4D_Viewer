%%
close all
clear all
clc

%%
x1 = linspace(0,1,21);
x2 = linspace(-1,1,41);
x3 = linspace(0,2,31);
x4 = linspace(-2,0,11);
[X1,X2,X3,X4] = ndgrid(x1,x2,x3,x4);
Y1 = X1 + sin(pi*X2) + X3.^2 + cos(pi*X4);


