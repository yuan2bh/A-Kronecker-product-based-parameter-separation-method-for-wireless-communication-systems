function [bq_est,bq_NorSquErr] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
% ��LchС�飬ÿС����Ե�һ��Ԫ�أ�����С���ȡƽ�����õ�PA�ķ��������ӹ���ʸ���������Ԫ��
% bq_est: bq�Ĺ���
% bq_NorSquErr: bq���Ƶ� ��һ����ƽ�������

if (length(chbq_est)~=Lch*Lbq)
    disp('errr! length(chbq_est)~=Lch*Lbq');
    return;
end
chbq_est_1=chbq_est(1 : Lbq );

bq_est =chbq_est_1 ./ chbq_est_1(1);

bq_err = bq_est-bq_CO;
bq_NorSquErr = (bq_err'*bq_err) / ( (bq_CO(2:end))'*(bq_CO(2:end)) );  %ȥ����һ��ϵ��1��Ӱ�죬��������ͬ

return;

