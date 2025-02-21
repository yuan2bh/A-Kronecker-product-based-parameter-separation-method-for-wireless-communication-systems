# A-Kronecker-product-based-parameter-separation-method-for-wireless-communication-systems

These MATLAB programs separate the nonlinear and linear parameters of nonlinear Hammerstein systems by a simple algorithm. They are the basis for our subsequent results. Although they were completed about 10 years ago and have some flaws, they are still valuable from today's point of view. 

The corresponding published Chinese papers are listed below: 
袁红林, 江立伟. 基于Kronecker积的无线通信系统参量分离方法[J].电讯技术, 2017, 57(10): 1099-1106. 
YUAN Honglin, JIANG Liwei. A Kronecker product based parameter separation method for wireless communication systems[J]. Telecommunication Engineering, 2017, 57(10): 1099-1106. (In Chinese) 

Abstract: For the problem of transmitter nonlinearity and wireless multipath fading channel in wireless communication systems, a combined estimation of the nonlinearity of power amplifier (PA) and impulse response of wireless channel with the received signal is proposed. Firstly, a convolution matrix is constructed based on the training symbols and model structure of the nonlinearity, then once least squares (LS) is used to estimate the coefficients of the model of the PA and actual transmitted symbols through the PA. Secondly, the estimation of the impulse response of the wireless channel is obtained with another LS. And the integration of the two estimations can be iterative or with direct summation and mean. Theoretical derivation and experimental results demonstrate that the method can be used in single-carrier or multi-carrier communication with oversampling technology, and can separate the nonlinearity of PA and impulse response of wireless channel effectively. The novel method can be used for the physical-layer radio frequency (ＲF) fingerprint authentication and reliable signal transmission.   
Key words: wireless channel; power amplifier nonlinearity; parameter separation; Kronecker product; radio frequency (ＲF) fingerprint identification; ＲF fingerprint authentication  

The used memory polynomial (MP) model of a nonlinear power amplifier (PA) with memory is described as
![image](https://github.com/user-attachments/assets/8df76466-2483-48fb-ba2b-ff53b1d5e586)  
(Equ. 2 in the published paper)  
where P and M are the nonlinear order and memory length respectively, and bp,m is the coefficients.

The main used functions:
function [s_cp] = obtain_train();
% The training sequence is based on the IEEE802.11a long training sequence samples, amplified and calibrated with an output power of 25dBm after amplification.

function [xn] = dn_WH_xn(dn,mem);
% Baseband input symbols after nonlinear Winner Hammerstein's (WH) PA are xn;
% dn: input;
% mem: 1 --- the PA with memory; 0 --- the PA without memory.

function [D,GMPmode] = P_MP_GMP_onlyOdd_D(dn,G,M,P);
% D: Nonlinear data matrix constructed based on the input dn of the GMP model;
% GMPmode: 'GMP', 'MP', 'P';  
% dn: input vector;
% G: delay depth; 
% M: memory depth; 
% P: nonlinear order. 

One of the entries is main_P.m. The run results correspond to Figure 7 in the published paper:  
![image](https://github.com/user-attachments/assets/096072f4-e6ac-4781-a4b5-f70e281afe7c)  
One result of main_P.m with MATLAB R2022a:  
![image](https://github.com/user-attachments/assets/9b78ba89-72a9-4a5c-a073-55b68104d7fa)  

The other entries is main_MP.m. The run results correspond to Figure 8 in the published paper:  
![image](https://github.com/user-attachments/assets/4e8c03a1-7308-4f80-bdc0-d9508f7475fa)  
One result of main_MP.m with MATLAB R2022a:  
![image](https://github.com/user-attachments/assets/3fa43b23-dafc-4bc9-a256-db38c15464a1)  


The results of the runs may not be identical, with possible reasons being the randomness of the noise, the lack of MontCarlo time or the version of Matlab etc. However, the trend of the results is the same.
