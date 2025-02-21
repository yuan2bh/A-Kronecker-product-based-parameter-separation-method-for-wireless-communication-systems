% run ok @ MATLAB R2022a 
% Honglin Yuan (Ԭ����)
% Nantong University, P.R.China (�й���ͨ��ѧ)
% 2025.2.20
%

% main.m
% GMP mode PA, separate nonlinearity coeff. and wireless channel
% for paper : һ�������豸�������������ŵ��ķ��뼼��
% 2015.4.3.
% yhl
% �������������PAϵ�����ܲ���Բ�����ʵPA���ݡ�
% 
% for MPģ�� of �����롷���� 

clear all;      clc;

% dn����ģ��PA������xn,���������ŵ�ch��Ϊy����AWGN������Ϊ�����ź�rn;����rn�����Ƶ��Ĺ�ʽ���Ƶõ������ŵ������������
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

% liu�ķ������1---transmitter1�� 2---transmitter2; 0---��������
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

%�����ŵ���λ������Ӧ;��MontCarlo���ȹ̶�
ch = [0.7; 0.3; 0.1]; 

Lch = length(ch);
    
dtLen=[1;2;4;8;16];  %ѵ�����и���

dtLenNum=length(dtLen);    
nmse_bq=zeros( length(snrVec),dtLenNum+1 );       %��һ��ΪSNR���ڶ���~��dtLenNum+1��Ϊ��Ӧ����ѵ�����е�nmse
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
        len_dn_gen= dtLen(mm);  %ѵ�����и���
        temp=obtain_train();        %25dBm;��MontCarlo���ȹ̶�
        dn = repmat(temp,len_dn_gen,1);
    end
    Ldn=length(dn);
    
    IBO=1.0
    dn=dn .* IBO;
    
    % dn����PA��
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

    % xn����channel��
    yn = conv(xn,ch);     %�����ź�ʸ��,xn��PA���
    Ly = length(yn);
    
        switch estMethod
            case 'Kalman'
                lenInit = 1*(Lbq*Lch);             % ��Kalman��ʼֵ���ݳ���

                len_iter=Ly;       % ����Kalman�������ݳ���

                bq_estMatMean = zeros(Lbq,len_iter);
                bq_NorSquErrVecMean = zeros(len_iter,1);
        end    
    
    % ��yn��һ��SNR��AWGN��������ÿһ��SNR��������Montcarlo���棻һ������ʵ�ֽ���һ�ι��ƣ�����ʵ�ֽ���mean�õ�һ��SNR��Ӧ�����
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
            % �������������Ӱ�첻��
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
                             chbq_estLS=pinv(D_D)*rn;     %α��, ������С���ˣ�����kroneker product : chbq
                             chbq_est=chbq_estLS;
                         otherwise
                             disp('LSmode is LS or SVD');
                             return;
                     end
                     
                     % ��chbq_est��LchС�飬ÿС����Ե�һ��Ԫ�أ�����С���ȡƽ�����õ�PA�ķ���������ʸ������bq_est�������Ĺ�һ����bq_err_sum
                     [bq_estMont(:,kk),bq_NorSquErrMont(kk)] = groupAveErr(chbq_est,Lch,Lbq,bq_CO);
                     
                     % ����bq�Ĺ��ƣ�����PA�������Ȼ������ŵ�
                     xn_est=D * bq_estMont(:,kk);    
                     xn_xn = X_X_gen(xn_est,Lch);   % ����PA���xn�ľ������
                     ch_est_temp=pinv(xn_xn) * rn;       % ���ŵ���LS����
                     ch_estMont(:,kk)=ch_est_temp;       % ÿ������ʵ�ֶ�����
                     
                     ch_err = ch_est_temp - ch;          % ���ŵ��������
                     ch_NorSquErrMont(kk) = (ch_err'*ch_err) / ( ch'* ch );   % ���ŵ����Ƶ����ģ�Ĺ�һ��(����ת�ã�)
                     
                 case 'Kalman'   
                     if (kk==1 & mm==1 & ii==1);   disp('estimate method is Kalman');           end
                     % ����Kalman��һ�����ݳ��ȡ�һ��SNR��ÿ������ʵ������������һ�ε������ѵ�����bq���Ƶ����ƽ����һ��ֵʸ������
                     bq_estMat=zeros(Lbq,len_iter);
                     bq_NorSquErrVec=zeros(len_iter,1);
                     [bq_estMat, bq_NorSquErrVec] = bq_Kalman(rn,D_D,Lbq,Lch,pNoise,lenInit,bq_CO);
                     
                     bq_estMatSum = bq_estMatSum + bq_estMat;
                     bq_NorSquErrVecSum = bq_NorSquErrVecSum + bq_NorSquErrVec;

                 otherwise
                     if (kk==1 & mm==1 & ii==1);   disp('the estimate method is beyond !');    return; end
             end
             
        end   % end MontCarlo of AWGN һ��SNR�µĲ�ͬnoiseʵ��
        
        % �Բ�ͬ���ݳ��ȡ���ͬSNR�µ�AWGN��montcarloʵ�ֺ�Ĵ���
        switch estMethod
            case  'LS'
                % ����LS���ƣ�������montcarloʵ�ֽ����ƽ��ֵ
                nmse_bq(ii,mm+1)= 10*log10( mean(bq_NorSquErrMont) );
                bq_est(ii,mm,:)=mean(bq_estMont,2);

                nmse_ch(ii,mm+1)= 10*log10( mean(ch_NorSquErrMont) );
                ch_est(ii,mm,:)=mean(ch_estMont,2);
                
            case 'Kalman'
                % Kalman������һ�����ݳ��ȡ�һ��SNR�²�ͬnoiseʵ�ֽ����ƽ��ֵ��Ȼ����bq_NorSquErrVecSum�е���Сֵor����ֵ�������Ӧ��bq����ʸ��
                bq_estMatMean = bq_estMatSum/montNum;
                bq_NorSquErrVecMean = bq_NorSquErrVecSum/montNum;
                
                % for display;��ÿһ�����ݳ��ȡ�ÿһ��SNR�µ�Montcarlo����ʵ�ֽ��ȡƽ������ʾ��
                figure;
                plot(10*log10(bq_NorSquErrVecMean));
                xlabel('SNR/dB'); ylabel('nmse of bq');
                title(strcat('symbols:',num2str(Ldn),'   SNR: ', num2str(snr),'dB'));
                
                [Y,I]=min(bq_NorSquErrVecMean);
                nmse_bq(ii,mm+1)= 10*log10( Y );
                bq_est(ii,mm,:)=bq_estMatMean(:,I);
            otherwise
                disp('wrong estimation method with �Բ�ͬ���ݳ��ȡ���ͬSNR�µ�AWGN��montcarloʵ�ֺ�Ĵ���'); 
        end
        
    end   %��ͬSNR
end   % ���ݳ��ȣ�headers��

% disp
bq_est(end,end,:)     %SNR=30;heasers=16ʱ��bq����ֵ

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




