clear;
clc;
close all;

%Script to read my c++ file
data_thing = importdata('output.txt');
x = data_thing(:,1)';
y = data_thing(:,2)';
plot(x,y,'*b')