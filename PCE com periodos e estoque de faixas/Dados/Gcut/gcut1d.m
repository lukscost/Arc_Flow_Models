%gcut1d

L =250; W = 250; re=0;

A=[167 184 90 0 0 
114 118 75 0 0
167 152 71 0 0
 83 140 13 0 0
 70  86 39 0 0
143 166 79 0 0
120 160 63 0 0
 66 148 93 0 0
 87 141 61 0 0
 69 165 85 0 0];

l=A(:,1);

w=A(:,2);

b=A(:,3);

% e = -10^(-10);
% 
% p = (max(l)*max(w)*100)/ (L*W)
% 
% diary on
% [AT, B, k, xb, x, FO, sr, E, arp] = GeracaoDeColunas(L,W,l,w, b, e,0);
% 
% %desenhaOti(AT,B,L,W,l,w,s,arp);
% diary off