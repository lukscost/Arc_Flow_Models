function [DE, U, Tt, gt, itTO, ft, o4] = MatMode(L,l,W,w,d,ri,mr)

% Essa funcao monta o modelo da Rita e do Valerio e resolve-o
%
% Entrada:
% L = tamanho do comprimento do objeto em estoque
% W = tamanho da largura do objeto em estoque
% l = vetor que contem a dimensï¿½o de cada comprimento dos itens
% w = vetor que contem a dimensï¿½o de cada largura dos itens
% d = vetor da demanda
% gerArc = 1 (geraÃ§Ã£o de colunas), diferente caso contrario

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
[w,I]=sort(w); l=l(I); d=d(I);
% quantidade de itens, sem rotação e com 


% roda a funcao Mat e MAtFa
[U, G, ga, ito] = Mat(W,w);
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
[lg,cg] = size(G);
[lu,cu] = size(U);
[lf,cf] = size(F);
[~,cd] = size(Dt);

% blocos nulos
o1 = zeros(lg,cf);
o2 = zeros(lu,cf -m);
o3 = zeros(lf,cg);

% cria a matriz da parte com igualdade do modelo
Tt = [G o1; U o2; o3 F];

% spy([Tt; DE],'ro',2), hold on
% spy(DE,'b+',2)

% se quiser vizualizar a matriz basta descomentar abaixo
% filename = 'MatrizIgualdade.txt';
% save('MatrizIgualdadeCompleta.mat','Tt')
% save MatrizIgualdadeCompleta.txt Tt -ascii

% para retirar as linhas nulas da matriz Tt
[liTt, lic] = size(Tt);
for i=1:liTt
       if isempty(find(Tt(i,:)))
       zerar(i,1) = i;
    end
    
end
[zerar,~,~] = find(zerar);  Tt(zerar,:) = [];

% save('MatrizIgualdadeSemLinhasNulas.mat','Tt')

% para saber a dimensao da matriz sem linhas nulas
[to4,rrr] = size(Tt);

% rhs da parte de igualdade
% if gerArc == 1
% o4 = zeros(to4,1);
% o4(1,1) = -1;
% o4(W+1,1) = 1;
% else
o4 = zeros(to4,1);
% end

% parte que cria as matriz da demandas
o5 = zeros(mr,cg +m);
DE = [o5 Dt];

% csvwrite('T0.csv',[full(Tt) full(o4); full(DE) d])

% se quiser vizualizar a matriz basta descomentar abaixo
% filename = 'matrizDemanda.xlsx';
% xlswrite(filename,DE,1,'A1')

end