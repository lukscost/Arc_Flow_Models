%gcut9d

L=1000; W=1000;

A=[310 426 91 0 0
673 468 61 0 0
426 463 40 0 0
325 498 95 0 0
555 540 18 0 0
292 455 30 0 0
343 341 56 0 0
362 491 95 0 0
305 688 51 0 0
459 607 53 0 0];

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