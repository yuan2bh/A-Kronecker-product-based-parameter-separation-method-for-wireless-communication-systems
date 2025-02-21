function [bq_estMat, bq_NorSquErrVec] = bq_Kalman(rn,D_D,Lbq,Lch,pNoise,lenInit,bq_CO);
% 采用 Kalman 估计 the kronecker product of wireless channel and nonlinearity coeff;
% bq_estMat: 迭代中的bq_est矢量矩阵，每一列为一个bq_est矢量
% bq_NorSquErrVec: 迭代中的bq估计误差平方归一化和矢量，对应bq_estMat中的bq_est矢量

% rn: 接收信号矢量
% D_D: 矩阵，由D根据信道长度右移并下降而成；是kronecker积的卷积矩阵
% Lbq: PA非线性系数总个数
% Lch：无线信道长度
% pNoise: power of noise
% lenInit: 基于LS求初始值的数据长度

% 使用前面数据基于LS获得kronecker product的初始值估计
chbq_init = pinv(D_D(1:lenInit,:)) * rn(1:lenInit);

% 使用剩下的数据基于Kalman进行迭代估计
% % r=rn(lenInit+1:end);
r=rn;

% obtain initial value accoring to page 257 of 何子述
temp1 = chbq_init - mean(chbq_init);
m = length(temp1)-1;
temp2 = corrmtx(temp1,m);
estCorr_init = temp2'*temp2;                  % 状态估计误差自相关矩阵的初始值 

L = length(r);
bq_estMat=zeros(Lbq,L);
bq_NorSquErrVec=zeros(L,1);
for n = 1: L
    % Kalman 迭代获得 channel 与 PA因子  kronecker product的估计
    if n==1                                   %step 1: 状态一步预测
        chbq_pred = chbq_init;         
    else
        chbq_pred=chbq_est;
    end
    alpha = r(n) - D_D(n,:) * chbq_pred;       % step 2: 根据观测计算新息
    if n==1
        preCorr = estCorr_init;                % step 3：一步预测误差自相关矩阵
    else
        preCorr = estCorr;
    end
    alphaCorr = D_D(n,:) * preCorr * D_D(n,:)'+ pNoise;      % step 4：新息过程自相关矩阵
    K = preCorr * D_D(n,:)'/ alphaCorr;                      % step 5: Kalman增益
    chbq_est = chbq_pred + K * alpha;                % step 6: state estimation
    estCorr = (eye(Lbq*Lch) - K * D_D(n,:)) * preCorr; % step 7：状态估计误差自相关矩阵
    
    %分组求bq等的估计,得到PA的非线性因子矢量估计bq_est及其误差的归一化和bq_NorSquErr
    [bq_estMat(:,n),bq_NorSquErrVec(n)] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
end


