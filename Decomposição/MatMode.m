function [A1eq, A1de, b1eq, A2eq, A2de, b2eq] = MatMode(L,l,W,w,ri,mr)

% Essa funcao monta o modelo da Rita e do Valerio e resolve-o
%
% Entrada:
% L = tamanho do comprimento do objeto em estoque
% W = tamanho da largura do objeto em estoque
% l = vetor que contem a dimens�o de cada comprimento dos itens
% w = vetor que contem a dimens�o de cada largura dos itens
% d = vetor da demanda
% gerArc = 1 (geração de colunas), diferente caso contrario

% Saida:
% F  = matriz das equacoes de balanceamento do fluxo de todos os s
% Dt = matriz para garantir que as demandas sejam satisfeitas
% gt = resultado da funcao Mat e MatFa (arcos)
% ites = resultado da funcao Mat e MatFa (itens)
% fs = resultado da funcao Mat e MatFa (faixa)


% PARA TESTES
%
% tipo = 'I';
%
% L=500;W=500;
% A=[198 205 91 0 0
% 179 155 71 0 0
% 364 236 94 0 0
% 272 147 74 0 0
% 352 145 32 0 0
% 343 245 63 0 0
% 132 174 60 0 0
% 164 250 27 0 0
% 282 356 58 0 0
% 342 151 75 0 0];
% l=A(:,1);
% w=A(:,2);
% d=A(:,3);
%
% gerArc = 0;
% W = 20; L = 30;
% w=[5;5;7;10;12];
% l=[7;10;12;8;10];
% d=[10;10;10;10;10];
%
% W = 11; L=11;
% w = [2;3;4;5]
% l = [2;3;4;5]
% d = [5;4;3;2]
%

% W = 15;
% L = 15;
%
% w=[2;3;4;5;7];
% l=[2;2;3;4;6];
% d=[10;10;10;10;10];

% para ordenar w em ordem crescente e a respectiva mudanca em l e d
% [w,I]=sort(w); l=l(I); d=d(I);

% roda a funcao Mat e MAtFa
[U, Go, ga, ito] = Mat(W,w);
[F,Dt,gi,ites,fs] = MatFa(L,l,w,ri,mr);

mm = length(unique(w));
% se quiser vizualizar a matriz basta descomentar abaixo
% filename = 'Dt.xlsx';
% xlswrite(filename,Dt,1,'A1')

% faz a juncao dos informacoes do grafo zero com os grafos s
gt = [ga gi];   itTO = [ito ites];
lito = length(ito);
itok = 8888*ones(1,lito);
oxe = lito - mm + 1;
itok(1)= 0;
itok(oxe:lito)=0;
% 9999 simboliza perda
[~,cito,~]=find(ito==9999);
itok(:,cito) = 0;
ft = [itok fs];

% quantidade de itens distintos e quantidade de itens
m = length(unique(w)); mt = length(w);

% para saber as dimensoes das matrizes e juntar os blocos
[lg,cg] = size(Go);
[lu,cu] = size(U);
[lf,cf] = size(F);
[ld,cd] = size(Dt);

zerar = [];
% matriz do segundo PLI
A1eq = Go;
[liTt, lic] = size(A1eq);
for i=1:liTt
    if isempty(find(A1eq(i,:)))
        zerar(i,1) = i;
    end
    
end
[zerar,~,~] = find(zerar);  A1eq(zerar,:) = [];
A1de = U(:,1:(cu-lu));

zerar = [];
% matriz do primeiro PLI
A2eq = F;
[liTt, lic] = size(A2eq);
for i=1:liTt
    if isempty(find(A2eq(i,:)))
        zerar(i,1) = i;
    end
    
end
[zerar,~,~] = find(zerar);  A2eq(zerar,:) = [];
A2de = [zeros(ld,lu) Dt];

% para saber a dimensao da matriz sem linhas nulas
[to4,rrr] = size(A1eq);
b1eq = zeros(to4,1);

[to4,rrr] = size(A2eq);
b2eq = zeros(to4,1);


end