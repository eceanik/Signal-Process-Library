function [no_of_sinusoid_mdl] = aksMusicMDL(y)

L = floor(length(y)/2);
X_beta = zeros(L,L);
for i=1:L
    
    X_beta(:,i) = y(i:i+L-1);
    
end;

Rs = (X_beta*X_beta')/L;




[~,d] = eig(Rs);
lamda = diag(d);

lamda = sort(lamda,'descend');
L = length(lamda);
N = 2*L-1;
% lamda = sort(lamda,'descend');
k_vec = 0:L-1;
mdl = zeros(1,length(k_vec));

for k  = 0:L-1
    num = 1;
    den = 0;
    for i = k+1:L
        num = num .* lamda(i)^(1/(L-k));
    end
    for i = k+1:L
        den = den + lamda(i);
    end
    liklihoodRatio = (num*(L-k)/den);
    logLiklihoodRatio = ((L-k)*N)*log(liklihoodRatio);
    mdl(k+1) = -logLiklihoodRatio + 0.5*k*(2*L - k)*log(N);
end
no_of_sinusoid_mdl = (find(mdl == min(mdl)));

end