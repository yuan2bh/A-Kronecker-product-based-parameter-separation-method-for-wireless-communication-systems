function [D_D] = D_D_gen(D,Ldn,Lbq,Lch);
% D根据Lch进行延迟与下沉，形成相对于   channel and 非线性因子的kronech积的 卷积矩阵: D_D
% rn: 接收信号矢量
% D: 矩阵，每一个行矢量是 基带符号（1个或多个） 根据GMP非线性模型结构构建
%      一个接收符号对应D的一行
% Ldn: 基带输入符号长度
% Lbq: PA非线性系数总个数
% Lch：无线信道长度
Lr=Ldn+Lch-1;
% Lr=length(rn);
D_D = zeros(Lr,Lbq*Lch);  
for i = 1:Lch
    D_D(i:i+Ldn-1,(i-1)*Lbq+1:i*Lbq) = D; 
end


