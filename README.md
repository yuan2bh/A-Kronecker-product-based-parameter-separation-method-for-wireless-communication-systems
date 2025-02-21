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

