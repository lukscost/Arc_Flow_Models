%gcut2d

L=250; W=250;
A=[ 120 133 61 0  0
135 186 75 0 0
 86  75 33 0 0
103  73 39 0 0
 66  85 20 0 0
135  97 79 0 0
 91 175 56 0 0
131 176 49 0 0
 71 176 34 0 0
153  72 96 0 0
 87 148 19 0 0
168 107 92 0 0
 18  90  8 0 0
140 109 79 0 0
132 159 13 0 0
152  93  2 0 0
135  68 78 0 0
121 158 40 0 0
 68  94 17 0 0
155  76 92 0 0];

l=A(:,1);

w=A(:,2);

b=A(:,3);

% e = -10^(-10);
% 
% p = (max(l)*max(w)*100)/ (L*W)
% 
% diary on
% [AT, B, k, xb, x, FO, sr, E, arp] = GeracaoDeColunas(L,W,l,w, b, e,1);
% 
% desenhaOti(AT,B,L,W,l,w,x,arp);
% diary off