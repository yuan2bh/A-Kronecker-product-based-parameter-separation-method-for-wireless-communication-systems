function [b_cal,nmse,GMPmode]=P_MP_GMP_onlyOdd_b_cal(G,M,P);
% P_MP_GMP_onlyOdd_b_cal.m
%  功能：产生非线性PA的P,MP与GMP模型的系数
%  dn:输入矢量
% P:非线性阶数, 奇数 or 偶数
% M:记忆深度
% G:delay深度
% 
%  yhl
%  2015.04.11
%

%训练序列头数
len_dn_gen = 30;     % bigger better  
% generate training sequence
temp = obtain_train();        %25dBm;可MontCarlo，先固定
dn = repmat(temp,len_dn_gen,1);

% P必须 >= 1
if (P<1)
    disp('P为>=1数');
    return;
end

if ( (G==0) && (M==0))
    GMPmode='P';
elseif ( (G==0) && (M~=0) )
        GMPmode='MP';
elseif ( (G~=0) && (M~=0) )
            GMPmode='GMP';
else
    disp('不考虑 G!=0 && M==0');
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

% 根据dn与GMP非线性模型参数G、M与P产生“输入非线性矩阵D”
[D,GMPmodeOut] = P_MP_GMP_onlyOdd_D(dn,G,M,P);

if(GMPmodeOut ~= GMPmode)
    disp('error! GMPmode~=GMPmodeOut');
    return;
end


% 根据是否有记忆(mem?=0)确定基带输入信号dn经过PA后的真实xn
xn = dn_WH_xn(dn,mem);


% 根据PA的实际输出xn估计PA在当前G、M、P参数(D)下的系数
b_est = pinv(D) * xn;       % 这里不能归一化，输出时归一化

% 根据PA的估计系数b_est得到在dn输入及当前G、M、P参数下的输出
xn_est = D * b_est;

% 在估计的模型系数下的输出误差性能，以此评判拟合模型是否准确
nmse = 10*log10( ( (xn_est-xn)' * (xn_est-xn)  ) ./ (xn' * xn) );

b_cal = b_est ./ b_est(1);
return;






