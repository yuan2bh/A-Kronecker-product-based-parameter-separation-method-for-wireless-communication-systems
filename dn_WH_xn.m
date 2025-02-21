function [xn] = dn_WH_xn(dn,mem);
% ����������ž���Winner Hammerstein (WH) PA��Ϊxn
% dn---H(z)---Amp+Fi---G(z)---xn
% �������ױ�д��
% dn: ��������
% mem: 0---�޼��䣻 1---�м���
% yhl
% 2015.4.11.

if mem==1
    % W_H �ĵ�һ���˲���H(Z)
    bH = [1; 0; 0.2];
    aH = [1;-0.1];
    temp1 = filter(bH,aH,dn);
elseif mem==0
    temp1 = dn;
end

% �����޼��������ϵͳTWTAģ��
temp1Amp = temp1 ./ ( 1 + 0.25 * temp1 .^ 2 );
temp1Fi = (0.26 * temp1 .^ 2) ./ ( 1 + 0.25 * temp1 .^ 2 );

temp2 = temp1Amp .* exp(i * temp1Fi);

if mem==1
    % W_H�ĵڶ����˲���G(Z)
    bG=[1; 0; -0.1];
    aG=[1; -0.2];
    xn=filter(bG, aG, temp2);
    % pCst(xn)
elseif mem==0
    xn = temp2;
end

return;


