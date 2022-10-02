function [Ve,Ite] = GSf(L,l,w,s,ri)

% Funcao que cria o grafo s (segundo estagio) do modelo da Rita e do Valerio de Carvalho,
% essa funcao usa somente os comprimentos do objeto e dos itens e a largura da faixa s
% criada no primeiro estagio, armazena os dados em uma matriz VE 
% ou Ite ambas (L x mve), onde mve eh o numero de itens que cabem na 
% faixa s, sendo que a se tem um arco em (i,i+l(j)) entao a matriz VE tem 
% na linha i+1 e na coluna j (do respectivo item) o valor de i+l(j) e Ite
% tem j ao invez de i+l(j)
%
% Entrada:
% L = tamanho do comprimento do objeto em estoque
% l = vetor que contem a dimensão de cada comprimento dos itens
% w = vetor que contem a dimensão de cada largura dos itens
% s = qual faixa eh
% 
% Saida:
% Ve = possui os arcos
% Ite = qual item aquele arco representa
%
%

% PARA TESTES

% clc
% clear all
% s=1;
% L = 30;
% w=[5;5;7;10;12];
% l=[7;10;12;8;10];

% 
% 
% s=3;
% L = 11;
% w = [2;3;4;5];
% l = [2;3;4;5];

% para ver quais itens podem ser usados na faixa
wu = unique(w);
[li,~,~] = find(w<=wu(s));

% vetor dos itens usados na faixa 's' e a quantidade de itens
ll=l(li);
mve = length(li);

% matriz dos arcos Ve e a matriz que informa qual item eh o arco
Ve = sparse(zeros(L,mve+1)); 
Ite = sparse(zeros(L,mve+1));

% para informar qual o maior arco que chegou em um vertice 
g = sparse(zeros(L,1));

% arco(s) do(s) item(s) inicias (da largura da faixa)
[li,~,~] = find(w==wu(s));
lp = ll(li);
Ve(1,li) = lp';
Ite(1,li) = ri(li,1);
% Ite(1,li) = li';
g(lp+1) = lp;

% para colocar todos os arcos possiveis depois dos iniciais
ss = length(li); % quantidade de itens da largura da faixa
for i = 1:ss 
for j = 1:mve

% se tiver somente um item da largura da faixa
if ss == 1
    it = lp(i) + ll(j); 
    if it <= L
        Ve(lp(i)+1,j) = it;
%         Ite(lp(i)+1,j)  = j;
        Ite(lp(i)+1,j)  = ri(j);
        g(it+1) = ll(j);
        if ll(j) > g(lp(i)+1) 
            g(lp(i)+1) = ll(j);
        end
    end

% se tiver mais que um item da largura da faixa
else
if ll(j) <= lp(i)
    it = lp(i) + ll(j); 
    if it <= L
        Ve(lp(i)+1,j) = it;
%         Ite(lp(i)+1,j) = j;
        Ite(lp(i)+1,j)  = ri(j);
        g(it+1) = ll(j);
        if ll(j) > g(lp(i)+1) 
            g(lp(i)+1) = ll(j);
        end
   end
end
end
end
end
% 
% full(Ve)
% full(Ite)

% EDITAR AQUI
% abaixo a ideia eh similar a apresentada para o grafo zero, 
% com exceção que soh considero como vertice aberto os vertices 
% apos o(s) inicia(s), porque nele nao se respeita a ordem  nao crescente
% dos arcos
[~,~,v] = find(Ve(2:L,:)); v=unique(v); tt = length(v);

v=sort(v);

if ~isempty(v)
vp = v(1);

for i= (vp+1):L

for j=1:mve
% HERE
    if ll(j) <= g(i)
        Q = i -1 + ll(j);
        if (Q <= L)
        Ve(i,j) = Q;
%         Ite(i,j) = j;
        Ite(i,j) = ri(j);
        g(Q+1) = ll(j);
            if ll(j) > g(Q+1) 
            g(Q+1) = ll(j);
            end
        end
    end 
end
end
end

% para colocar as perdas depois do vertice min(w)+1
% full(Ve)
% p = min(Ve(1,1:mve));
p = ll(s);

% mostra onde houve os arcos de perdas
indp = p:(L-1);


% coloca a informacao em Ve e em Ite (perda sera representada por 9999)
Ve(indp+1,mve+1) = indp+1;
Ite(indp+1,mve+1) = 9999;

% full(Ve)
% full(Ite)
% end