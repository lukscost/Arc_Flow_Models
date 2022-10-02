% Costa, L. L. S. (2022), Extens�es do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matem�tica Estat�stica e Computa��o Cient�fica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-est�gios e com aparo.
% Programa para problema com ou sem rota��o

% Resumo da adaptacao do modelo cem rota��o:
% - Rotaciono soh os itens retangulares e atualizo o conjunto dos itens;
% - Crio um vetor 'ri' que informa qual item eh.
% Ex: (3,4) -> ri_1 = 3, (4,3) -> ri_2 = 3 (mesmo indice 'ri'),
% ou seja, itens rotacionados tem o mesmo valor em 'ri';
% - Na criacao da matriz dos grafos do segundo est�gio 
% cria-se uma matriz correspondente 'Ite' que informa os indices dos itens;
% - Na cria��o da 'Ite' uso o 'ri' para infomar o indice real do item;
% - A matriz 'Ite' � usada na cria��o da matriz da demanda para 
% informar o indice do item;
% - As demais coisas s�o iguais.

clc
clear all
close all
warning off 
diary('Output')

% com rota��o: hotrool = 1 
% sem rota��o: caso contr�rio
hotrool = 0;

if hotrool == 1
disp('Problema COM rota��o de itens')
else
disp('Problema SEM rota��o de itens')
end

% caminho do cplex
opl

dork = [];
amax = [];
funT = [];
peT =  [];
cys = {};

mkdir('CSV')

for po = 1:12 

% dados
por = ['Gcut\gcut',num2str(po),'d.m'];
run(por);
[w,I]=sort(w); l = l(I); d = d(I); 

disp(['Executando problema: ',num2str(po),' ...'])

% quantidade de itens, sem rota��o e com 
mwf = length(w);

if hotrool == 1
% com rota��o
% para nao duplicar itens quadrados
Id = find((w-l)==0);
rota = (1:mwf)';
rota(Id) = [];

% para rotacionar os itens que nao sao quadrados
ri = [(1:mwf)'; rota];
wr = [w; l(rota)];
lr = [l; w(rota)];
dr = [d; d(rota)];
else
% sem rota��o
ri = (1:mwf)';
wr = w;
lr = l;
dr = d;    
end

[wr,I]=sort(wr); lr=lr(I); ri = ri(I); dr = dr(I); 
mr = length(wr);

% problema, W, L, itens, demanda 
dri = [po W L mwf mr min(wr) max(wr)  min(lr) max(lr)  min(d) max(d)];
dork = [dork; dri];

% Matriz do modelo
try
tic
disp('Iniciando montagem da matriz do modelo ...')
[DE, U, Tt, gt, itTO, ft, o4] = MatMode(L,lr,W,wr,ri,mwf);
tia = toc; 
Erro = 0;

catch
disp('Erro: N�o foi poss�vel criar a matriz do modelo')
Erro = 1;
tia = inf;
end

if Erro == 1
    continue
end

[mt,nt] = size(Tt);
[mdel,ndel] = size(DE);
[m1,n1] = size([Tt; DE]);
[mmax,nmax] = Ordem(W,L,lr,wr,d);
% (%) of the elements in the matrix are zeros.
dens = (numel([Tt; DE])-nnz([Tt; DE]))*100/numel([Tt; DE]);
amax0 =[po m1 mmax n1 nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')
% vetor de custos do PI
c = zeros(nt,1);
c(1) = 1;

% Para gerar os limitantes
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = 8*60; %8 min

[x,fo,u,y] = cplexlp(c,-DE,-d,Tt,o4,LB,UB,[],opt);
if isempty(fo)
    y;
    disp('Erro: Cplex n�o resolveu o problema relaxado.')
    continue
end
I = find(c);
x(I);
y;
fo;

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');

[xu,foi,u,yi] = cplexmilp(c,-DE,-d,Tt,o4,[],[],[],LB,UB,type,[],opt);
cys = [cys; [num2str(po), ') ', yi.cplexstatusstring]]; 
if isempty(foi)
    yi;
    disp('Erro: Cplex n�o resolveu o problema inteiro.')
    continue
end
% quantidade de objetos usados
obj = xu(I);
yi;
foi;

funp = [po fo foi tia yi.time];
funT = [funT; funp];

% avalia algumas porcentagens
if ~isempty(xu)
    r = DE*xu;
    tok = (L.*W)'*xu(I);
    pe = tok - ((r')*(l.*w));
    perda  = (pe*100)/tok;
    exce = (sum(r) - sum(d))*100;
    exc = exce/sum(d);
    foi = fo;
else    
    disp('ERROR')
    foi  = 0;
    perda = 10000;
    exc = 10000;
% end
% y
end

% informa��es das perdas
pep = [po perda exc  (sum(r) - sum(d))];
peT = [peT; pep];
disp(' ' )
disp('Itens, demanda-i, produ��o-i, i=1,2,3')
legdi = ['Demanda e produ��o; Colunas: 1) Itens; '...
    ' 2) Demanda; 3) Produ��o;'];
tit = 'Itens. & $d_{i}$ & $pro_{i}$';
ProtoP = [(1:mwf)' d  r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(po);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)
end


disp('========================================================' )
disp('RECORTAR AQUI')

disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = ['Informa��es do problema; Colunas: 1) Problema; '...
    ' 2) Valor da largura do objeto, W;' ...
    ' 3) Valor do comprimento do objeto, L;' ...
    ' 4) Quantidade de itens sem rota��o;'...
    ' 5) Quantidade de itens com rota��o;'...
    ' 6) Intervalo da largura dos itens (w);'...
    ' 7) Intervalo do comprimento dos itens (l);'...
    ' 8) Intervalo da demanda de itens (d).'];
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 0 1 0 1 0 1 0])
csvwrite('CSV/WLwld.csv',dork)


disp(' ' )
disp(' ' )
disp('Dimens�o da matriz:') 
legdi = ['Dimens�o da matriz do modelo; Colunas: 1) Problema; '...
    ' 2) N�mero de restri��es (nr);' ...
    ' 3) Limitante no n�mero de restri��es (lnr);' ...
    ' 4) N�mero de vari�veis; (nv);' ...
    ' 5) Limitante no n�mero de vari�veis (lnv);'...
    ' 6) $\%$ de elementos nulos (nz).'];
tit = 'Prob. & nr & lnr & nv & lnv & nz';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amax,legdi,tit,[ ],1,[])
csvwrite('CSV/mnModelo.csv',amax)

disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = 'Prob. & f.o.r & f.o.i & T. mod. & T. Cplex ';
legdi = ['Dados da solu��o; Colunas: 1) Problema; '...
    ' 2) Valor da fun��o-objetivo do problema relaxado (f.o.r.);' ...
    ' 3) Valor da fun��o-objetivo do problema inteiro (f.o.i);' ...
    ' 4) Tempo para constru��o da matriz do modelo;'...
    ' 5) Tempo para resolu��o do modelo.'];
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)

disp(' ' )
disp(' ' )
disp('desperdicio da �rea total \% , exec produ \%, exec em itens') 
legdi = ['An�lise da produ��o; Colunas: 1) Problema; '...
    ' 2) Porcentagem de desperd�cio da �rea total dos objetos;' ...
    ' 3) Porcentagem do excesso de produ��o da quantidade total de itens produzidos;' ...
    ' 4) Quantidade de itens produzidos a mais.'];
tit = 'Prob. & $(\%$ - desperd�cio$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)

disp(cys)

diary off