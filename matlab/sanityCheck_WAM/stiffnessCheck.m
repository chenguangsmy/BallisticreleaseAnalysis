clear all
close all
clc

%syms ks kr fp xs xr x1 x2
%syms ks xs 
%syms kr xr
syms x2
fp = 15;
x1 = 0;
%x2 = -0.02734;
ks = 320;
kr = 300;
xs = 0.05;
xr = -0.0533;
%xr = 0.482;

% Sanity check with force before release
%f1 = kr*(xr-x1)
 
% Solve for stiffness and x0
% [ks xs] = solve([0  == ks*(xs-x1)+kr*(xr-x1),...
%                  fp == ks*(xs-x2)+kr*(xr-x2)],[ks xs]);
%[kr, xr] = solve([0  == ks*(xs-x1)+kr*(xr-x1),...
%                 fp == ks*(xs-x2)+kr*(xr-x2)],[kr xr]);
%[x1 x2] = solve([0  == ks*(xs-x1)+kr*(xr-x1),...
%                 fp == ks*(xs-x2)+kr*(xr-x2)],[x1 x2]);
x2 = solve(fp == ks*(xs-x2)+kr*(xr-x2), x2)
double(x2)
%double(kr)
%double(xr)
             
