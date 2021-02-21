function [f,Px] = aksMusic(y,Fs,L,fRange,p)


% f_lower = fRange(1);
% diff = fRange(2);
% f_upper = fRange(3);

f_lower = fRange.f_lower;
f_upper = fRange.f_upper;
diff = fRange.diff;

f = f_lower:diff:f_upper;
w = f*2*pi/Fs;
A = zeros(length(w),L);

% X_beta = zeros(L,L);
% 
% for i=1:L
%     
%     X_beta(:,i) = y(i:i+L-1);
%     
% end
[~,Rs] = corrmtx(y,L-1,'autocorrelation');


% Rs = (X_beta*X_beta')/L;
% p = aksMusicMDL(Rs);
% p = 4;


for i=1:L-p
    
        A(:,i) = exp(-1*1i*w*(i-1))';
%     A(:,i) = cos(w*(i-1))';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v,d] = eig(Rs);
lamda = diag(d);
[~,ind]  = sort(lamda,'descend');

j = p+1:L;
noiseEigVec = v(:,ind(j));


B = A*noiseEigVec;
B = abs(B);
B = B.^2;

Px = sum(B,2);
Px = 1./Px;
end