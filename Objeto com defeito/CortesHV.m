% codigo que cria os novos objetos após retirar os objetos
function [Oh,Ov,des] = CortesHV(W,L,wd,ld,XY) 
% clc
% clear all
% close all
% warning off
% 
% diary('teste')
% % dados:
% % Tamanho do objeto
% W = 200;
% L = 250;
% % coordenda x,y (começando de baixo na esquerda do objeto)
% XY = [99 99];
% % tamanho do defeito
% wd = 1;
% ld = 2;

% objetos obtidos por cortes na horizontal
Wh1 = XY(1,2);
Wh2 = W - XY(1,2) - wd;
Lh = L;
Whk = [Wh1; Wh2];
Lhk = [L; L];

% objetos obtidos por cortes na horizontal
Wv = W;
Lv1 = XY(1,1);
Lv2 = L - XY(1,1) - ld;
Wvk = [W; W];
Lvk = [Lv1; Lv2];

% objetos finais
Oh = [Whk Lhk];
Ov = [Wvk Lvk];

% desperdicio: area da faixa horizontal e vertical
des = [wd*L  ld*W ];













