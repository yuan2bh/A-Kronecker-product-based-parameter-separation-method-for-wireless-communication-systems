function [D_D] = D_D_gen(D,Ldn,Lbq,Lch);
% D����Lch�����ӳ����³����γ������   channel and ���������ӵ�kronech���� �������: D_D
% rn: �����ź�ʸ��
% D: ����ÿһ����ʸ���� �������ţ�1�������� ����GMP������ģ�ͽṹ����
%      һ�����շ��Ŷ�ӦD��һ��
% Ldn: ����������ų���
% Lbq: PA������ϵ���ܸ���
% Lch�������ŵ�����
Lr=Ldn+Lch-1;
% Lr=length(rn);
D_D = zeros(Lr,Lbq*Lch);  
for i = 1:Lch
    D_D(i:i+Ldn-1,(i-1)*Lbq+1:i*Lbq) = D; 
end


