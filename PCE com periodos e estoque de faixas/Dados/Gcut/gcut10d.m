%gcut10d

L=1000; W=1000;

A=[730 300 68 0 0
269 717 47 0 0
463 369 69 0 0
642 464  2 0 0
453 329 78 0 0
455 667 82 0 0
506 482 42 0 0
560 362 52 0 0
483 260 13 0 0
693 345 32 0 0
381 510  8 0 0
456 586 43 0 0
457 453 12 0 0
707 658 12 0 0
639 650 93 0 0
691 359 12 0 0
434 700 21 0 0
576 291 39 0 0
728 739  6 0 0
545 742 99 0 0];

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