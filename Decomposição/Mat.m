function [U, G, ga, ito] = Mat(W,w)

% Essa funcao recebe o valor da largura do objeto e as larguras dos itens e
% cria a matriz G e U do modelo, a matriz G eh a matriz de equacao de
% balanceamento do fluxo e a matriz U serve para garantir que as faixas
% contruidas do primeiro estagio serao usadas
% 
% Entrada:
% W = tamanho da largura do objeto em estoque
% w = vetor que contem a dimens�o de cada largura dos itens
% 
% Saida:
% G = matriz das equacoes de balanceamento do fluxo
% U = matriz para garantir que as faixas sejam usadas no segundo estagio
% ga = matriz que informa que na coluna j da matriz G, eh referente ao
%      arco ga(2,j) - ga(1,j), isto eh, ga' contem os arcos do grafo G0
% ito = vetor que informa na sua posicao j qual item eh o arco da coluna j
%
%

% PARA TESTES
%
% clc
% clear all
% W = 11;
% w = [2;3;4;5];
% 
%  W = 20;
%  w=[5;5;7;10;12];

% roda a funcao que cria a matriz dos arcos
[Ve, Ite ] = GO(W,w);
Ve = sparse(Ve);
Ite = sparse(Ite);

% filename = 'Grafo0gcut2.xlsx';
% xlswrite(filename,Ve,1,'A1')
% xlswrite(filename,Ite,2,'A1')

% filename = 'GrafoS.xlsx';
% xlswrite(filename,Ve,pag-1,'A1')
% xlswrite(filename,Ite,pag,'A1')

% filename = 'Grafo0C5.xlsx';
% Ve = xlsread(filename,1);
% Ite = xlsread(filename,2);



nv = nnz(Ve); % quantidade de arcos

% inicia a matriz G e a coluna da variavel z (fluxo do grafo zero)
G = sparse(zeros(W+1,nv+1));
G(1,1)=1;
G(W+1,1)=-1;


m = length(unique(w)); % quantidade de itens de largura distintas

% inicia a matriz ga, ito e U
ga = sparse(zeros(2,nv+2+m));
ito = sparse(zeros(1,nv+2+m));
U = sparse(zeros(m,nv+1+m));

% criacao das colunas de zs que garantem que as faixas serao usadas
for i=1:m
    U(i,nv+1+i) = -1;    
end

% preenchimento das colunas de G e U
n = 2;
for i=1:W
    for j=1:m
      
% percorro a matriz Ve e quando encontro um arco (uma entrada nao nula)
% crio uma coluna para esse arco na matriz G
        if Ve(i,j)~=0

% se o arco nao eh de perda informo na matriz U que aquela coluna eh de 
% uma faixa que serah usada       
         if Ite(i,j) ~= 9999
         G(i,n) = -1;
         G(Ve(i,j)+1,n) = 1;

         ga(1,n) = i;
         ga(2,n) = Ve(i,j)+1;
         
         U( Ite(i,j),n) = 1;
         
         ito(1,n) = Ite(i,j); % faixa 0 serah representada como item 8888

         else
         G(i,n) = -1;
         G(i+1,n) = 1;
         ga(1,n) = i;
         ga(2,n) = i;
         ito(1,n) = 9999;
         end
                  
         n = n + 1;       
        end
    
    end
end

% arcos de perda
Dei = eye(W-min(w));
er = zeros(1,W-min(w));
De1 = [er; Dei];
De2 = [-Dei; er];
Def = De1 + De2;

for i=1:min(w)
   Def = [er; Def];
end



% +2 pq tem o fluxo
Iny = nv -(W-min(w))+2;
G(:,Iny:end) = Def;

Iny2 = Iny + (W-min(w)) -1;
ga(:,Iny:Iny2) = [(min(w)+1:W); (min(w)+2:W+1)];
ito(:,Iny:Iny2) = 9999;

ga(:,end) = [];
ito(:,end) = [];

Ing = find(ga(1,:));


% full([ga' ito'])
% full([(ga(:,Ing)0)' ito(:,Ing)'])
% Grafos(W, [(ga(:,Ing)-1)' ito(:,Ing)'], [], [], [])

% full(ga)
% full(G)

% se quiser salvar a matriz no excel para vizualizar eh soh descomentar
% filename = 'EqBalGrafo0.xlsx';
% xlswrite(filename,G,1,'A1')
% xlswrite(filename,U,2,'A1')
%
% U
% full(G)
% ga
% ito

% end