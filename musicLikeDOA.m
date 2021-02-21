function [f,Px,B,Rs] = musicLikeDOA(y,Fs,L,fRange)

T = 1/Fs;

X_beta = zeros(L,L);
f_lower = fRange.f_lower;
f_upper = fRange.f_upper;
diff = fRange.diff;

w = window(@chebwin,L);
for i=1:L
    
    X_beta(:,i) = y(i:i+L-1).*w';
%     X_beta(:,i) = y(i:i+L-1);
    
end
% Rs = (X_beta*X_beta')/(L*L);
% Rs = (X_beta*X_beta')/L;
Rs = (X_beta'*X_beta)/L;      % As per Thesis


f = f_lower:diff:f_upper;
w = f*2*pi*T;
A = zeros(length(w),L);
 
 for i=1:L
     
     A(:,i) = exp(-1*1i*w*(i-1))';
%      A(:,i) = cos(w*(i-1))';
     
 end
% A = transpose(A);             % As per Thesis


B = A*Rs*A';
% B = A'*Rs*A;                  % As per Thesis
B = abs(B);
Px = diag(B)/L;
% Px = Px/L^3;









end