clear all;
% close all;
clc;
fs=100;
f=1; % Hz
t=[0:1/fs:100]';
y=sin(2*pi*f*t);
fv=[0.1:0.1:25];
[TFR,timeVec,freqVec]=morletTFR(y,fv,fs,3);
[X,Y]=meshgrid(timeVec, freqVec);
figure;
Z = zeros([length(X) 1]);
imagesc(timeVec,freqVec,TFR)
colorbar

