%gcut5d

L=500;W=500;

A=[198 205 91 0 0
179 155 71 0 0
364 236 94 0 0
272 147 74 0 0
352 145 32 0 0
343 245 63 0 0
132 174 60 0 0
164 250 27 0 0
282 356 58 0 0
342 151 75 0 0];

l=A(:,1);

w=A(:,2);

b=A(:,3);

% e = -10^(-10);
% p = (max(l)*max(w)*100)/ (L*W)
% 
% diary on
% [AT, B, k, xb, x, FO, sr, E, arp] = GeracaoDeColunas(L,W,l,w, b, e,1);
% 
% %desenhaOti(AT,B,L,W,l,w,s,arp);
% diary off