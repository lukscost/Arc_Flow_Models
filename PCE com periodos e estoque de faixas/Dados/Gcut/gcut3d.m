%gcut3d

L=250; W=250; 

A=[ 66  80 33 0 0
164 107 75 0 0
 64 184 94 0 0
121  86 65 0 0
163 135  2 0 0
 85  98 78 0 0
 81 102 49 0 0
103 186  4 0 0
152 106  7 0 0
176 139  6 0 0
111 118  5 0 0
 69 169 20 0 0
146 133 26 0 0
112 112 55 0 0
133 160 50 0 0
 63 129 64 0 0
163 152 25 0 0
110 155 80 0 0
 96 136 90 0 0
 92 142 84 0 0
 84 143 95 0 0
119 133 33 0 0
 71  71 64 0 0
146  84  5 0 0
 93  86 68 0 0
 89  86 30 0 0
101 146 100 0 0
172  73 33 0  0
 73 169 62 0  0
 99  99 87 0  0];

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
% %desenhaOti(AT,B,L,W,l,w,s,arp);
% diary off