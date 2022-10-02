function [m,n] = Ordem(W,L,l,w,wr)

wu = unique(w);
twu =length(wu);
% linhas da matriz da demanda
tl =length(wr);

% para ver quais itens podem ser usados na faixa
wu = unique(w);

% 1 estagio
varz0 = 1;
qa1e = ((W+1)^2-(W+1))/2;
% 2 estagio
varzs = length(wu);
qa2e = 0;

for i=1:twu
    [li,~,~] = find(w==wu(i));
    ll=l(li);
    lt = min(ll);
    qa2e = qa2e + (((L-lt+1)^2-(L-lt+1))/2) +1; 
end

% QUANTIDADE MAXIMA DE COLUNAS (VARIAVEIS)
n = varz0 + qa1e + varzs + qa2e;

% QUANTIDADE MAXIMA DE LINHAS (RESTRICOES)
m = (W+1) + varzs + varzs*(L+1) + tl;
end