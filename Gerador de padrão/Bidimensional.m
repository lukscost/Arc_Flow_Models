% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios,
% Modelo gerador de padrão de corte IRRESTRITO e que permite ou não rotação

clc
clear all
close all
warning off 

% script para salvar o caminho das funções do CPLEX
opl

% se rot = 1, permite rotação,
% outro valor é caso contrário.
rot = 1;

if rot == 1
    disp('Modelo gerador de padrão de corte irrestrito e COM rotação de itens')
    diary('TestesCOM')
else
    disp('Modelo gerador de padrão de corte irrestrito e SEM rotação de itens')
    diary('TestesSEM')
end
mkdir('CSV')

% matrizes dos dados
dork = []; amax = []; funT = [];
peT =  []; cys = {};


for po = 1:1
    
disp(['Executando problema: ',num2str(po),' ...'])
    
% load dos dados
por = ['Dados/gcut',num2str(po),'d.mat'];
load(por);
[w,I]=sort(w); l = l(I); d = d(I); 
[w l];

% quantidade de itens, sem rotação e com 
mwf = length(w);

if rot == 1
% com rotação
% para nao duplicar itens quadrados
Id = find((w-l)==0);
rota = (1:mwf)';
rota(Id) = [];
ri = [(1:mwf)'; rota];
wr = [w; l(rota)];
lr = [l; w(rota)];
dr = [d; d(rota)];
else
% sem rotação
ri = (1:mwf)';
wr = w;
lr = l;
dr = d;
end

[wr,I]=sort(wr); lr=lr(I); ri = ri(I); dr = dr(I); 
% [wr lr ri dr];

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
disp('Erro: Não foi possível criar a matriz do modelo')
Erro = 1;
tia = inf;
end

if Erro == 1
    continue
end

% Para gerar os limitantes na dimensão da matriz
[mt,nt] = size(Tt);
[m1,n1] = size([Tt; DE]);
[mmax,nmax] = Ordem(W,L,lr,wr,w);
dens = (numel([Tt; DE])-nnz([Tt; DE]))*100/numel([Tt; DE]);
amax0 =[po mt mmax (nt-1) nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

% vetor de custos 
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
    disp('Erro: Cplex não resolveu o problema relaxado.')
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
    disp('Erro: Cplex não resolveu o problema inteiro.')
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

funp = [po -fo -foi tia yi.time perda sum(r) ];
funT = [funT; funp];


% [W L]
% [w l d r] % r é o padrão de corte
disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = ['Demanda e produção; Colunas: 1) Itens; '...
    ' 2) Demanda; 3) Produção;'];
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
legdi = 'Informações do problema';
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 0 1 0 1 0 1 0])
csvwrite('CSV/WLwld.csv',dork)


disp(' ' )
disp(' ' )
disp('Dimensão da matriz:') 
legdi = 'Dimensão da matriz do modelo';
tit = 'Prob. & nr & lnr & nv & lnv & nz';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amax,legdi,tit,[ ],1,[])
csvwrite('CSV/mnModelo.csv',amax)


disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = 'Prob. & f.o.r & f.o.i & T. mod. & T. Cplex  & perda & itens';
legdi = 'Dados da solução';   
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)

disp(cys)

diary off
