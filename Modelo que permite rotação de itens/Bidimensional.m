% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios,
% Modelo de fluxo em arcos que permite a rotação de 90º

% Resumo da adaptacao do modelo com rotação em relação ao tradicional:
% - Rotaciono soh os itens retangulares e atualizo o conjunto dos itens;
% - Crio um vetor 'ri' que informa qual item eh.
% Ex: (3,4) -> ri_1 = 3, (4,3) -> ri_2 = 3 (mesmo indice 'ri'),
% ou seja, itens rotacionados tem o mesmo valor em 'ri';
% - Na criacao da matriz dos grafos do segundo estágio 
% cria-se uma matriz correspondente 'Ite' que informa os indices dos itens;
% - Na criação da 'Ite' uso o 'ri' para infomar o indice real do item;
% - A matriz 'Ite' é usada na criação da matriz da demanda para 
% informar o indice do item;
% - As demais coisas são iguais.

clc
clear all
close all
warning off 
diary('Testes')

mkdir('CSV')

% script para salvar o caminho das funções do CPLEX
opl

% matrizes dos dados
dork = []; amax = []; funT = [];
peT =  []; cys = {};

for po = 1:1 
    
disp(['Executando problema: ',num2str(po),' ...'])

% load dos dados
por = ['Gcut\gcut',num2str(po),'d.mat'];
load(por);
[w,I]=sort(w); l = l(I); d = d(I); 
% [w l d];

% quantidade de itens, sem rotação e com 
mwf = length(w);

% para nao duplicar itens quadrados
Id = find((w-l)==0);
rota = (1:mwf)';
rota(Id) = [];

% para rotacionar os itens que nao sao quadrados
ri = [(1:mwf)'; rota];
wr = [w; l(rota)];
lr = [l; w(rota)];
dr = [d; d(rota)];
%[wr lr ri dr];
mr = length(wr);

[wr,I]=sort(wr); lr=lr(I); ri = ri(I); dr = dr(I); 
% [wr lr ri dr];

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
disp('Erro: Não foi possível criar a matriz do modelo')
Erro = 1;
tia = inf;
end
if Erro == 1
    continue
end

% Para gerar os limitantes na dimensão da matriz
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
tipe = 'I'*ones(nt,1);
type = char(tipe');
[xu,foi,u,yi] = cplexmilp(c,-DE,-d,Tt,o4,[],[],[],LB,UB,type,[],opt);
cys = [cys; [num2str(po), ') ', yi.cplexstatusstring]]; 
if isempty(foi)
    yi;
    disp('Erro: Cplex não resolveu o problema inteiro.')
    continue
end
% quantidade de objetos usados
obj = xu(I);
yi;
foi;

% função objetivo e tempos computacionais
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

% informações das perdas
pep = [po perda exc  (sum(r) - sum(d))];
peT = [peT; pep];

disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = ['Demanda e produção; Colunas: 1) Itens; '...
    ' 2) Demanda; 3) Produção;'];
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
tit = 'Prob. & f.o.r & f.o.i & T. mod. & T. Cplex ';
legdi = 'Dados da solução';
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)


disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = 'Análise da produção';
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)

disp(cys)