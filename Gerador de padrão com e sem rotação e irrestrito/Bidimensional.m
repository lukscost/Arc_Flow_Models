% Costa, L. L. S. (2022), Extens�es do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matem�tica Estat�stica e Computa��o Cient�fica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-est�gios e com aparo.
% Programa para o modelo gerador de padr�o de corte 
% com e sem rota��o para o caso irrestrito

clc
clear all
close all
warning off 
diary('Output')
 
% com rota��o: hotrool = 1 
% sem rota��o: caso contr�rio
hotrool = 1;

if hotrool == 1
disp('Modelo gerador de padr�o de corte irrestrito e COM rota��o de itens')
else
disp('Modelo gerador de padr�o de corte irrestrito e SEM rota��o de itens')
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
% po = 1;

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
dri = [po W L mwf mr min(w) max(w)  min(l) max(l)  min(d) max(d)];
dork = [dork; dri];

% Matriz do modelo
try
tic
disp('Iniciando montagem da matriz do modelo ...')
[DE, U, Tt, gt, itTO, ft, o4, idz, id1es, idzs , id2es] = MatMode(L,lr,W,wr,ri,mwf);
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
[m1,n1] = size([Tt; DE]);
[mmax,nmax] = Ordem(W,L,lr,wr,w);
% (%) of the elements in the matrix are zeros.
dens = (numel([Tt; DE])-nnz([Tt; DE]))*100/numel([Tt; DE]);
amax0 =[po mt mmax (nt-1) nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

% vetor de custos do PI
c0 = full(itTO');
c0(c0==9999) = inf;
c0(c0==8888) = inf;

for uy = 1:mwf
    Ind = find(itTO==uy);
    c0(Ind,1) = l(uy)*w(uy);
end

c = zeros(nt,1);
c(id2es,1) = c0(id2es,1);
c(c==inf) = 0;

% Para gerar os limitantes
LB = sparse(zeros(nt-1,1));
UB = inf*ones(nt-1,1);

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = 8*60; %8 min

% atualiza as matrizes
Ma = Tt(:,2:end);
rhs = -Tt(:,1);
c = c(2:end,1);

[x,fo,u,y] = cplexlp(-c,[],[],Ma,rhs,LB,UB,[],opt);
% cys = [cyz; y.cplexstatusstring]; 
if isempty(fo)
    y;
    disp('Erro: Cplex n�o resolveu o problema relaxado.')
    continue
end
% quantidade de objetos usados
I = find(c);
x(I);
y;
fo;

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt-1,1);
type = char(tipe');

[xu,foi,u,yi] = cplexmilp(-c,[],[],Ma,rhs,[],[],[],LB,UB,type,[],opt);
cys = [cys; [num2str(po), ') ', yi.cplexstatusstring]]; 
if isempty(foi)
    yi;
    disp('Erro: Cplex n�o resolveu o problema inteiro.')
    continue
end
% quantidade de objetos usados
obj = 1;
yi;
foi;

% avalia algumas porcentagens
if ~isempty(xu)
    r = DE*[0; xu];
    tok = L*W;
    pe = tok - ((r')*(l.*w));
    perda  = (pe*100)/tok;
    foi = fo;
else    
    disp('ERROR')
    foi  = 0;
    perda = 10000;
    exc = 10000;
% end
% y
end

funp = [po -fo -foi tia yi.time perda sum(r)];
funT = [funT; funp];
disp(' ' )
disp('Itens, demanda-i, produ��o-i, i=1,2,3')
legdi = ['Demanda e produ��o; Colunas: 1) Itens; '...
    ' 2) Demanda; 3) Produ��o;'];
tit = 'Itens. & $d_{i}$ & $pro_{i}$';
ProtoP = [(1:mwf)' w l  r];
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
tit = 'Prob. & f.o.r & f.o.i & T. mod. & T. Cplex & Perda & itens ';
legdi = ['Dados da solu��o; Colunas: 1) Problema; '...
    ' 2) Valor da fun��o-objetivo do problema relaxado (f.o.r.);' ...
    ' 3) Valor da fun��o-objetivo do problema inteiro (f.o.i);' ...
    ' 4) Tempo para constru��o da matriz do modelo;'...
    ' 5) Tempo para resolu��o do modelo'...
    ' 6) Desperd�cio do padr�o;'...
    ' 7) Quantidades de itens no padr�o de corte;'];
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)

disp(cys)

diary off