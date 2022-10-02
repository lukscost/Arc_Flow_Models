function [F,Dt,gi,ites,fs] = MatFa(L,l,w,ri,mr)

% Essa funcao recebe o valor do comprimento do objeto e dos itens e
% cria a matriz F e Dt do modelo, a matriz F eh a matriz blocada que 
% contem as equacoes de balamceamento do fluxo e cada bloco eh referemte a
% um grafo s, a matriz Dt serve para garantir que as demandas sejam
% satisfeitas
%
% Entrada:
% L = tamanho do comprimento do objeto em estoque
% l = vetor que contem a dimensão de cada comprimento dos itens
% w = vetor que contem a dimensão de cada largura dos itens
% 
% Saida:
% F  = matriz das equacoes de balanceamento do fluxo de todos os s
% Dt = matriz para garantir que as demandas sejam satisfeitas
% ga = mesmo da funcao Mat
% ites = vetor que informa na sua posicao j qual item eh o arco da coluna j
% fs = vetor que informa qual faixa s a coluna j pertence
%

% PARA TESTES
% 
% clc
% clear all
% L = 11;
% w = [2;3;4;5];
% l = [2;3;4;5];
% 
% L = 30;
% w=[5;5;7;10;12];
% l=[7;10;12;8;10];


% quantidade de itens distintos e quantidade de itens
m = length(unique(w)); mt = length(w);

% iniciacao das matrizes usadas
F = sparse(zeros(m*L+m*1,m+1)); Dt = []; gi = []; fs = []; ites = []; na=m;

% inicio do loop
for s=1:m

% crio a matriz dos arcos da faixa s    
[Ve,Ite] = GSf(L,l,w,s,ri);
Ve = sparse(Ve);
Ite = sparse(Ite);

% pag = 2*s;
% filename = 'GrafoSgcut12.xlsx';
% xlswrite(filename,Ve,pag-1,'A1')
% xlswrite(filename,Ite,pag,'A1')

% C5
% filename = 'GrafoSgcut2.xlsx';
% sheet = pag-1;
% Ve = xlsread(filename,sheet);
% Ite = xlsread(filename,pag);

% para prencher as entradas de zs na matriz
fim = s*L+s*1;  ini = fim - L;
F(ini,s)=+1;    F(fim,s)=-1;

% numero de elementos nao nulos de Ve e dimensao
nv = nnz(Ve);   [~,m] = size(Ve);

% parte que inicia a matriz G (que contem a equacao de balanceamento de fluxo 
% do grafo s - faixa s), ga (arcos), Ds (que informa quais itens s pode
% satisfazer na demanda), fa (informa a faixa), itea (informa o item)

G = sparse(zeros(L+1,nv));  ga = sparse(zeros(2,nv));    Ds = sparse(zeros(mr,nv));   
fa = sparse(zeros(1,nv));  itea = sparse(zeros(1,nv));

% loops
n=1;
for i=1:L
    for j=1:m
% se o valor de Ve for diferente de zero, ou seja, se possui um arco no 
% grafo entao se prenche uma coluna do bloco e a informacao nas outras
% matrizes
        if Ve(i,j)~=0
         G(i,n) = -1;
         G(Ve(i,j)+1,n) = +1;
         
         % se nao for perda
         if Ite(i,j) ~= 9999
         ga(1,n) = i;
         ga(2,n) = Ve(i,j)+1;
         fa(1,n) = s;
         itea(1,n) = Ite(i,j);
         dsi = Ite(i,j);
         Ds(dsi,n) = 1;
         
         %se for perda
         else
         ga(1,n) = i;
         ga(2,n) = i+1;
	     fa(1,n) = s;
	     itea(1,n) = Ite(i,j);
         end
         
         % para ir para a proxima coluna do bloco
         n = n + 1;       
        end
    
    end
end


% valores para preencher a matriz que contem os blocos nas possicoes
% corretas
nc = na+1;
na = na + nv;
nf = nc:na;

nl = s*L + s;
ni = nl - L;
nn = ni:nl;

% size(G)
% Gg = full(G);
% save  Gg;

[~,Ind] = sort(ga(1,:));

Ds = Ds(:,Ind);
ga = ga(:,Ind);
fa = fa(:,Ind);
itea = itea(:,Ind);
G = G(:,Ind);

% colocando os blocos nas matrizes 
Dt = [Dt Ds];   gi = [gi ga];   fs = [fs fa];   ites = [ites itea];
F(nn,nf) = G;

% full(ga)

Ing = find(ga(1,:));
% full([(ga(:,Ing)-1)' itea(:,Ing)'])
% Grafos(L, [(ga(:,Ing)-1)' itea(:,Ing)'], [], [], [])

% pause

% se quiser salvar em excel basta descomentar abaixo
% filename = 'EqBalGrafoS.xlsx';
% xlswrite(filename,G,s,'A1')

end

% filename = 'matriz.xlsx';
% xlswrite(filename,Dt,m+3,'A1')
% xlswrite(filename,F,m+4,'A1')

% se quiser vizualizar a matriz basta descomentar abaixo
% F
% Dt
% gi
% ites
% fs

end