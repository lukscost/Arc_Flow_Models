function [m,n] = Ordem(W,L,l,w)

wu = unique(w);
twu =length(wu);
tl =length(l);

% para ver quais itens podem ser usados na faixa
wu = unique(w);



% QUANTIDADE MAXIMA DE COLUNAS (VARIAVEIS)
varz0 = 1;
varzs = length(wu);

qa1e = ((W+1)^2-(W+1))/2;
qa2e = 0;

for i=1:twu
    [li,~,~] = find(w==wu(i));
    ll=l(li);
    lt = min(ll);
    qa2e = qa2e + ((L-lt+1)^2-(L-lt+1))/2 +1; 
end

n = varz0 + qa1e + varzs + qa2e + tl;


% QUANTIDADE MAXIMA DE LINHAS (RESTRICOES)
m = (W+1) + varzs + varzs*(L+1) + tl;
end