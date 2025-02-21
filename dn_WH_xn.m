function [xn] = dn_WH_xn(dn,mem);
% 基带输入符号经过Winner Hammerstein (WH) PA后为xn
% dn---H(z)---Amp+Fi---G(z)---xn
% 根据文献编写。
% dn: 输入数据
% mem: 0---无记忆； 1---有记忆
% yhl
% 2015.4.11.

if mem==1
    % W_H 的第一个滤波器H(Z)
    bH = [1; 0; 0.2];
    aH = [1;-0.1];
    temp1 = filter(bH,aH,dn);
elseif mem==0
    temp1 = dn;
end

% 经过无记忆非线性系统TWTA模型
temp1Amp = temp1 ./ ( 1 + 0.25 * temp1 .^ 2 );
temp1Fi = (0.26 * temp1 .^ 2) ./ ( 1 + 0.25 * temp1 .^ 2 );

temp2 = temp1Amp .* exp(i * temp1Fi);

if mem==1
    % W_H的第二个滤波器G(Z)
    bG=[1; 0; -0.1];
    aG=[1; -0.2];
    xn=filter(bG, aG, temp2);
    % pCst(xn)
elseif mem==0
    xn = temp2;
end

return;


