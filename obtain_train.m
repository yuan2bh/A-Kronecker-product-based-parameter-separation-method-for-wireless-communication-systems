function [s_cp] = obtain_train()
% 
% 该训练序列是根据IEEE802.11a的long训练序列样本，放大经过功放后输出功率为25dBm标定而成；
% 程序为：LTrain160_64_25dBm_IN.m
% yhl
% 2015.1.28
%
load('LTrain160_64_25dBm_IN.mat','LTrain160_64_25dBm_IN'); 
s_cp=LTrain160_64_25dBm_IN;