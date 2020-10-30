function [Ve, Ite ] = GO(W,w)

% Funcao que cria o grafo zero (primeiro estagio) do modelo da Rita e do Valerio de Carvalho,
% essa funcao usa somente a largura do objeto e a largura dos itens, ela
% armazena os dados em uma matriz VE ou Ite (W x wu), onde wu é o tamanho dos elemnetos
% unicos de w, sendo que a se tem um arco em (i,i+w(j)) entao a matriz VE tem 
% na linha i+1 e na coluna j (do respectivo item) o valor de i+w(j) e Ite
% tem j ao invez de i+w(j)
%
% Entrada:
% W = tamanho da largura do objeto em estoque
% w = vetor que contem a dimensão de cada largura dos itens
% 
% Saida:
% Ve = possui os arcos
% Ite = qual item aquele arco representa
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

w=unique(w); % w sem repeticoes 

m = length(w); % m eh quantos itens nao repetidos se tem em w

% Ve = sparse(zeros(W+1,m+1)); % indica pra onde vão os vertíces
Ve = zeros(W+1,m+1); % indica pra onde vão os vertíces
Ite = zeros(W+1,m+1);

% coloca os arcos iniciais (0,w(j)) que ficara em (1,w(j))
% porque no Matlab indice comeca do 1
Ve(1,:) = [w' 0];
Ite(1,:) = [1:m 0];

g = zeros(W+1,1); % indicarah qual o indice maior arco que o vertice recebeu
g(w+1) = 1:m; % arcos iniciais

ta = zeros(W+1,1); % indicarah qual o maior arco que o vertice recebeu
ta(w+1) = w; % arcos iniciais

tic
for i = (w(1)+1):(W+1)
    in = g(i);
    for j = 1:in
        Q = i + w(j)-1;
        if (Q <= W) % se colocar o item nao excede a dimensao da largura 

        Ve(i,j) = Q; % entao coloco o item
        Ite(i,j) = j; % e guardo a posicao
        
        if w(j) <= g(Q) || g(Q)==0  % se item eh menor que o ultimo arco posto
            if g(Q) < j % atualizo o tamanho do arco que aquele vertice 
               g(Q) = j; % recebeu se for necessario
            end
        end
        end

    end
end



% para colocar as perdas depois do vertice min(w)+1
p = min(w);

% mostra onde houve os arcos de perdas
indp = p:(W-1);

% coloca a informacao em Ve e em Ite (perda sera representada por 9999)
Ve(indp+1,m+1) = indp+1;
Ite(indp+1,m+1) = 9999;

Ve(end,:) = [];
Ite(end,:) = [];

% a = toc
% [(0:W-1)' Ve]
% [(0:W-1)' Ite]
% Ve
% Ite

% end
