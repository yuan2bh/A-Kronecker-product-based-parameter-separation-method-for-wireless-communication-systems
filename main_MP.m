% run ok @ MATLAB R2022a 
% Honglin Yuan (袁红林)
% Nantong University, P.R.China (中国南通大学)
% 2025.2.20
%

% main.m
% GMP mode PA, separate nonlinearity coeff. and wireless channel
% for paper : 一种无线设备非线性与无线信道的分离技术
% 2015.4.3.
% yhl
% 由于随机产生的PA系数性能差，所以采用真实PA数据。
% 
% for MP模型 of 《分离》论文 

clear all;      clc;

% dn经过模型PA的数据xn,经过无线信道ch后为y，与AWGN噪声后为接收信号rn;采用rn根据推导的公式估计得到无线信道与非线性因子
format long;

estMethod = 'LS'
% estMethod = 'Kalman'

% for debug program
DEBUG= 0

% LSmode= 'SVD'
LSmode= 'LS'

%Montecarlo Simulation times
if DEBUG    montNum = 1
else        montNum = 100;
end
disp(['monNum is: ',num2str(montNum)]);

if DEBUG     snrVec=[50]
else         snrVec=[5:5:30]
end    

% % % % 'GMP', 'MP', 'P'
% GMPmode= 'P' 
GMPmode= 'MP' 
% GMPmode= 'GMP' 

% liu的发射机，1---transmitter1； 2---transmitter2; 0---不是刘的
LIU = 0

switch GMPmode
    case 'P'
        G=0
        M=0
        P=7
        if (LIU==1 || LIU==2)
            P=7
        end
    case 'MP'
        G=0
        M=3
        P=7
    case 'GMP'
        G=2
        M=4
        P=7
    otherwise
        disp('mode is : P,MP,GMP!');
        return;
end

%无线信道单位脉冲响应;可MontCarlo，先固定
ch = [0.7; 0.3; 0.1]; 

Lch = length(ch);
    
dtLen=[1;2;4;8;16];  %训练序列个数

dtLenNum=length(dtLen);    
nmse_bq=zeros( length(snrVec),dtLenNum+1 );       %第一列为SNR，第二列~第dtLenNum+1列为相应长度训练序列的nmse
nmse_ch=zeros( length(snrVec),dtLenNum+1 );  

switch GMPmode
    case 'P'
        bq_est=zeros(length(snrVec),dtLenNum,(P+1)/2);
    case 'MP'
        bq_est=zeros(length(snrVec),dtLenNum, ( (P+1)/2)*(M+1) );
    case 'GMP'
        bq_est=zeros(length(snrVec),dtLenNum, ((P+1)/2)*(M+1)  + ((P+1)/2-1)*(M+1)*G );
    otherwise 
        disp('wrong GMP mode in bq_est');
        return;       
end

ch_est=zeros(length(snrVec),dtLenNum,Lch);

for mm=1:dtLenNum
    if DEBUG
        len_dn_gen= 5
        rng(34)
        dn = rand(len_dn_gen,1)+i*rand(len_dn_gen,1);
    else
        len_dn_gen= dtLen(mm);  %训练序列个数
        temp=obtain_train();        %25dBm;可MontCarlo，先固定
        dn = repmat(temp,len_dn_gen,1);
    end
    Ldn=length(dn);
    
    IBO=1.0
    dn=dn .* IBO;
    
    % dn经过PA后
    [D,GMPmode_D] = P_MP_GMP_onlyOdd_D(dn,G,M,P);
    [bq_CO,nmseCal,GMPmode_bq_CO]=P_MP_GMP_onlyOdd_b_cal(G,M,P);  %calibration

    if (GMPmode=='P')
        switch LIU
            case 1
                bq_CO=[1;-0.0735-0.0114*i; -0.0986+0.059*i; -0.0547-0.0055*i]; %liu, Ttr1 @ P mode
            case 2
                bq_CO=[1; -0.091+0.158*i; 0.2503+0.0286*i; 0.0155+0.0025*i];   %liu, Ttr2 @ P mode
            otherwise
                disp('P mode PA is not from LIU, generated with WH mode');
        end
    end
    
    Lbq = length(bq_CO);    
    [D_D] = D_D_gen(D,Ldn,Lbq,Lch);
    
    if (mm==1);   
        rank_D_D=rank(D_D)  
        rank_D=rank(D)
     end
    
    xn=D*bq_CO;     %PA output
    
    chbq = kron(ch,bq_CO);     
    Lchbq=length(chbq);

    % xn经过channel后
    yn = conv(xn,ch);     %接收信号矢量,xn是PA输出
    Ly = length(yn);
    
        switch estMethod
            case 'Kalman'
                lenInit = 1*(Lbq*Lch);             % 求Kalman初始值数据长度

                len_iter=Ly;       % 用于Kalman迭代数据长度

                bq_estMatMean = zeros(Lbq,len_iter);
                bq_NorSquErrVecMean = zeros(len_iter,1);
        end    
    
    % 对yn加一定SNR的AWGN噪声，对每一个SNR噪声进行Montcarlo仿真；一个噪声实现进行一次估计，所有实现进行mean得到一个SNR对应结果。
    for ii = 1:length(snrVec)       % SNR(dB)
        snr = snrVec(ii);
        nmse_bq(ii,1)=snr;
        nmse_ch(ii,1)=snr;

        switch estMethod
            case 'Kalman'
                bq_estMatSum=zeros(Lbq,len_iter);
                bq_NorSquErrVecSum=zeros(len_iter,1);
        end
        
        bq_NorSquErrMont=zeros(montNum,1);
        bq_estMont=zeros(Lbq,montNum);
        ch_estMont=zeros(Lch,montNum);
        ch_NorSquErrMont=zeros(montNum,1);
                     
        for kk = 1:montNum   %MontCarlo of AWGN 
            % 随机噪声对性能影响不大
            if DEBUG            rng('default'); 
            else
                rng('default'); 
                rng('shuffle');     % use current time as seed of random numbers;
            end

             rn=awgn(yn,snr,'measured');
             
             noise = rn-yn;
             pNoise=var(noise);

             switch estMethod
                 case 'LS'
                     if (kk==1 & mm==1 & ii==1);   disp('estimate method is LS');   end

                     switch LSmode
                         case 'SVD'
                             D_D_squ = D_D' * D_D;
                             [U,S,V]=svd(D_D_squ);
                             sizeV=size(V);
                             chbq_estSVD = zeros(Lchbq,1);
                             door=length(D_D_squ) * norm(D_D_squ) *eps(class(D_D_squ));
                             for iii=1:sizeV(2)
                                if ( S(iii,iii)> door ) 
                                    chbq_estSVD = chbq_estSVD + V(:,iii)' * D_D' * rn * (1/S(iii,iii)) * V(:,iii);  
                                end
                             end
                            chbq_est=chbq_estSVD;
                         case 'LS'
                             chbq_estLS=pinv(D_D)*rn;     %伪逆, 采用最小二乘，估计kroneker product : chbq
                             chbq_est=chbq_estLS;
                         otherwise
                             disp('LSmode is LS or SVD');
                             return;
                     end
                     
                     % 把chbq_est分Lch小组，每小组除以第一个元素，所有小组加取平均，得到PA的非线性因子矢量估计bq_est及其误差的归一化和bq_err_sum
                     [bq_estMont(:,kk),bq_NorSquErrMont(kk)] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
                     
                     % 根据bq的估计，计算PA的输出；然后估计信道
                     xn_est=D * bq_estMont(:,kk);    
                     xn_xn = X_X_gen(xn_est,Lch);   % 产生PA输出xn的卷积矩阵
                     ch_est_temp=pinv(xn_xn) * rn;       % 求信道的LS估计
                     ch_estMont(:,kk)=ch_est_temp;       % 每次噪声实现都保存
                     
                     ch_err = ch_est_temp - ch;          % 求信道估计误差
                     ch_NorSquErrMont(kk) = (ch_err'*ch_err) / ( ch'* ch );   % 求信道估计的误差模的归一化(共轭转置！)
                     
                 case 'Kalman'   
                     if (kk==1 & mm==1 & ii==1);   disp('estimate method is Kalman');           end
                     % 对于Kalman，一个数据长度、一个SNR的每个噪声实现样本，都有一次迭代，把迭代中bq估计的误差平方归一化值矢量保存
                     bq_estMat=zeros(Lbq,len_iter);
                     bq_NorSquErrVec=zeros(len_iter,1);
                     [bq_estMat, bq_NorSquErrVec] = bq_Kalman(rn,D_D,Lbq,Lch,pNoise,lenInit,bq_CO);
                     
                     bq_estMatSum = bq_estMatSum + bq_estMat;
                     bq_NorSquErrVecSum = bq_NorSquErrVecSum + bq_NorSquErrVec;

                 otherwise
                     if (kk==1 & mm==1 & ii==1);   disp('the estimate method is beyond !');    return; end
             end
             
        end   % end MontCarlo of AWGN 一个SNR下的不同noise实现
        
        % 对不同数据长度、不同SNR下的AWGN的montcarlo实现后的处理
        switch estMethod
            case  'LS'
                % 对于LS估计，求噪声montcarlo实现结果的平均值
                nmse_bq(ii,mm+1)= 10*log10( mean(bq_NorSquErrMont) );
                bq_est(ii,mm,:)=mean(bq_estMont,2);

                nmse_ch(ii,mm+1)= 10*log10( mean(ch_NorSquErrMont) );
                ch_est(ii,mm,:)=mean(ch_estMont,2);
                
            case 'Kalman'
                % Kalman，先求一定数据长度、一定SNR下不同noise实现结果的平均值；然后找bq_NorSquErrVecSum中的最小值or收敛值，及其对应的bq估计矢量
                bq_estMatMean = bq_estMatSum/montNum;
                bq_NorSquErrVecMean = bq_NorSquErrVecSum/montNum;
                
                % for display;对每一个数据长度、每一个SNR下的Montcarlo噪声实现结果取平均，显示：
                figure;
                plot(10*log10(bq_NorSquErrVecMean));
                xlabel('SNR/dB'); ylabel('nmse of bq');
                title(strcat('symbols:',num2str(Ldn),'   SNR: ', num2str(snr),'dB'));
                
                [Y,I]=min(bq_NorSquErrVecMean);
                nmse_bq(ii,mm+1)= 10*log10( Y );
                bq_est(ii,mm,:)=bq_estMatMean(:,I);
            otherwise
                disp('wrong estimation method with 对不同数据长度、不同SNR下的AWGN的montcarlo实现后的处理'); 
        end
        
    end   %不同SNR
end   % 数据长度，headers数

% disp
bq_est(end,end,:)     %SNR=30;heasers=16时的bq估计值

ch_est(end,end,:)

figure;
subplot(211);
plot(nmse_bq(:,1),nmse_bq(:,2:end),'-.*');
xlabel('SNR(dB)');
ylabel('nmse - PA');
grid on;
legend('headers:1','headers:2','headers:4','headers:8','headers: 16');
% title('PA coefficients');
title('with memory(MP)');

subplot(212);
plot(nmse_ch(:,1),nmse_ch(:,2:end),'-.*');
xlabel('SNR(dB)');
ylabel('nmse - wireless channel');
grid on;
legend('headers:1','headers:2','headers:4','headers:8','headers: 16');
% title('wireless channel');

return;




