function [X_X] = X_X_gen(X,Lch);
% X����Lch�����ӳ����³����γ������   channel �� �������: X_X
% Lch�������ŵ�����
Lx=length(X);
Lr=Lx+Lch-1;
X_X = zeros(Lr,Lch);  
for i = 1:Lch
    X_X(i:i+Lx-1,i) = X; 
end


