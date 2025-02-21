function [b_cal,nmse,GMPmode]=P_MP_GMP_onlyOdd_b_cal(G,M,P);
% P_MP_GMP_onlyOdd_b_cal.m
%  ���ܣ�����������PA��P,MP��GMPģ�͵�ϵ��
%  dn:����ʸ��
% P:�����Խ���, ���� or ż��
% M:�������
% G:delay���
% 
%  yhl
%  2015.04.11
%

%ѵ������ͷ��
len_dn_gen = 30;     % bigger better  
% generate training sequence
temp = obtain_train();        %25dBm;��MontCarlo���ȹ̶�
dn = repmat(temp,len_dn_gen,1);

% P���� >= 1
if (P<1)
    disp('PΪ>=1��');
    return;
end

if ( (G==0) && (M==0))
    GMPmode='P';
elseif ( (G==0) && (M~=0) )
        GMPmode='MP';
elseif ( (G~=0) && (M~=0) )
            GMPmode='GMP';
else
    disp('������ G!=0 && M==0');
    return;
end

switch GMPmode
    case 'P'
        mem=0
    case 'MP'
        mem=1
    case 'GMP'
        mem=1
    otherwise
        disp('GMPmose is :G, MP,or GMP');
        return;
end

% ����dn��GMP������ģ�Ͳ���G��M��P��������������Ծ���D��
[D,GMPmodeOut] = P_MP_GMP_onlyOdd_D(dn,G,M,P);

if(GMPmodeOut ~= GMPmode)
    disp('error! GMPmode~=GMPmodeOut');
    return;
end


% �����Ƿ��м���(mem?=0)ȷ�����������ź�dn����PA�����ʵxn
xn = dn_WH_xn(dn,mem);


% ����PA��ʵ�����xn����PA�ڵ�ǰG��M��P����(D)�µ�ϵ��
b_est = pinv(D) * xn;       % ���ﲻ�ܹ�һ�������ʱ��һ��

% ����PA�Ĺ���ϵ��b_est�õ���dn���뼰��ǰG��M��P�����µ����
xn_est = D * b_est;

% �ڹ��Ƶ�ģ��ϵ���µ����������ܣ��Դ��������ģ���Ƿ�׼ȷ
nmse = 10*log10( ( (xn_est-xn)' * (xn_est-xn)  ) ./ (xn' * xn) );

b_cal = b_est ./ b_est(1);
return;






