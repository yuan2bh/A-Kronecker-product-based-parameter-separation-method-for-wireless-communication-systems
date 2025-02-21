function [bq_estMat, bq_NorSquErrVec] = bq_Kalman(rn,D_D,Lbq,Lch,pNoise,lenInit,bq_CO);
% ���� Kalman ���� the kronecker product of wireless channel and nonlinearity coeff;
% bq_estMat: �����е�bq_estʸ������ÿһ��Ϊһ��bq_estʸ��
% bq_NorSquErrVec: �����е�bq�������ƽ����һ����ʸ������Ӧbq_estMat�е�bq_estʸ��

% rn: �����ź�ʸ��
% D_D: ������D�����ŵ��������Ʋ��½����ɣ���kronecker���ľ������
% Lbq: PA������ϵ���ܸ���
% Lch�������ŵ�����
% pNoise: power of noise
% lenInit: ����LS���ʼֵ�����ݳ���

% ʹ��ǰ�����ݻ���LS���kronecker product�ĳ�ʼֵ����
chbq_init = pinv(D_D(1:lenInit,:)) * rn(1:lenInit);

% ʹ��ʣ�µ����ݻ���Kalman���е�������
% % r=rn(lenInit+1:end);
r=rn;

% obtain initial value accoring to page 257 of ������
temp1 = chbq_init - mean(chbq_init);
m = length(temp1)-1;
temp2 = corrmtx(temp1,m);
estCorr_init = temp2'*temp2;                  % ״̬�����������ؾ���ĳ�ʼֵ 

L = length(r);
bq_estMat=zeros(Lbq,L);
bq_NorSquErrVec=zeros(L,1);
for n = 1: L
    % Kalman ������� channel �� PA����  kronecker product�Ĺ���
    if n==1                                   %step 1: ״̬һ��Ԥ��
        chbq_pred = chbq_init;         
    else
        chbq_pred=chbq_est;
    end
    alpha = r(n) - D_D(n,:) * chbq_pred;       % step 2: ���ݹ۲������Ϣ
    if n==1
        preCorr = estCorr_init;                % step 3��һ��Ԥ���������ؾ���
    else
        preCorr = estCorr;
    end
    alphaCorr = D_D(n,:) * preCorr * D_D(n,:)'+ pNoise;      % step 4����Ϣ��������ؾ���
    K = preCorr * D_D(n,:)'/ alphaCorr;                      % step 5: Kalman����
    chbq_est = chbq_pred + K * alpha;                % step 6: state estimation
    estCorr = (eye(Lbq*Lch) - K * D_D(n,:)) * preCorr; % step 7��״̬�����������ؾ���
    
    %������bq�ȵĹ���,�õ�PA�ķ���������ʸ������bq_est�������Ĺ�һ����bq_NorSquErr
    [bq_estMat(:,n),bq_NorSquErrVec(n)] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
end


