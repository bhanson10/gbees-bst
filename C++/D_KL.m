% Kullback-Liebler divergence
% as a test vs the C++ code
%and also allows doing it for the analytical diffusion case solved in matlab

function [H,Hm] = D_KL(P,Q)

% assume / assert same sizes

% normalise
P(isnan(P))=0;
P = P/sum(sum(P));
Q(isnan(Q))=0;
Q = Q/sum(sum(Q));

Hm = real(P.*log(P./Q));
Hm(find(isnan(Hm)))=0; % 0*log(0) = 0
H=sum(sum(Hm));
H=H/log(2);
% has miniscule imag part due to error
