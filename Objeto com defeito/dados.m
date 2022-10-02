% PARA TESTES: Objetos com defeitos

% dados:
% Tamanho do objeto
W0 = 100;
L0 = 200;

% jt = 1;
orw = jt;

XYk =[ 100    50
   100    40
   100    60
   125    20
   125    71
    30    30
    80    40
    80    40];

% coordenda x,y (começando de baixo na esquerda do objeto)
XY = XYk(orw,:); %[comprimento L; altura W]

lwk =[5 4;
     5 4 ;
     5 6 ;
    7 10 ;
    7 8  ;
   10 10 ;
   30 18 ;
   38 18 ];

% tamanho do defeito
ld = lwk(orw,1);
wd = lwk(orw,2);

[Oh,Ov,des] = CortesHV(W0,L0,wd,ld,XY);

% WL = [Oh; Ov; W0 L0]
WL = [Oh; Ov]

w = [30    26    20    35    22]';
l = [40    68    50    60    45]';
v = [10    12     8    18     9]';
d = [10; 10; 10; 10; 10];

% [w,I]=sort(w); l=l(I); d = d(I); v = v(I);

% [(1:5)' w l v]

