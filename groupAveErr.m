function [bq_est,bq_NorSquErr] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
% 分Lch小组，每小组除以第一个元素，所有小组加取平均，得到PA的非线性因子估计矢量及其误差元素
% bq_est: bq的估计
% bq_NorSquErr: bq估计的 归一化、平方、误差

if (length(chbq_est)~=Lch*Lbq)
    disp('errr! length(chbq_est)~=Lch*Lbq');
    return;
end
chbq_est_1=chbq_est(1 : Lbq );

bq_est =chbq_est_1 ./ chbq_est_1(1);

bq_err = bq_est-bq_CO;
bq_NorSquErr = (bq_err'*bq_err) / ( (bq_CO(2:end))'*(bq_CO(2:end)) );  %去除第一个系数1的影响，与文献相同

return;

