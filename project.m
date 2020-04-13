clear all; close all; clc;

sessions{1} = load('data/Mouse12-120806_awakedata.mat');
sessions{2} = load('data/Mouse28-140313_awakedata.mat');

for i = 1:length(sessions)
   data = sessions{i};
   
end