clc
clear all
close all

%%

worldmap('world')
c = load('coast');
plotm(c.lat,c.long,'k')
