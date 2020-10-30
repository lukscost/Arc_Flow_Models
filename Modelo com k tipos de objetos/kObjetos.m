% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios,
% com aparo e múltiplos tipos de objetos

clc
clear all
close all
warning off 
format long

% script para salvar o caminho das funções do CPLEX
opl

diary('TestesGcut')
csvdir = 'CSVGcutXdv';
mkdir(csvdir)

% inicia as matrizes de dados
fou = [];  datis = [];
amak = []; dork = [];
cys = {}; peT = [];
ort = [];

disp('Ax >= d')

for jt= 1:1
gei = [];

disp(['Executando problema: ',num2str(jt),' ...'])

% load dos dados
datjt = ['DataGcut/gcut',num2str(jt),'dv.mat'];
load(datjt)
k = length(W);
m = length(w);

disp(['Dimensão: ',num2str(k),'x',num2str(m)])

%  organização dos dados
[w,I]=sort(w); l=l(I); d = d(I);
[mo,no] = size(W);

% inicia as matrizes do modelo
PatK = [];
DemK = [];
rhs = [];
ck = [];
amai = [0 0];

disp('Iniciando montagem da matriz do modelo ...')
try
tic;
for i=1:mo

% função responsável pela criação das matrizes de cada tipo de objeto
[DE, U, Tt, gt, itTO, ft, I, o4, tew, cd, Dt, idz, id1es, idzs , id2es] = MatMode(L(i,1),l,W(i,1),w,d);

[mt,nt] = size(PatK);
[mk,nk] = size(Tt);

if mt == 0 % PARA GERAR O BLOCO K=1
    PatK = Tt;
    DemK = DE;
    rhs = o4;
    
    % vetor de custos
    c = zeros(nk,1);
    c(1) = L(i,1)*W(i,1);
    ck = c;
else
    PatK = [PatK zeros(mt,nk); zeros(mk,nt) Tt];
    DemK = [DemK  DE];
    rhs = [rhs; o4];
    
    % vetor de custos
    c = zeros(nk,1);
    c(1) = L(i,1)*W(i,1);
    ck = [ck; c];
end

[mmax,nmax] = Ordem(W(i,1),L(i,1),l,w);
amai0 = [mmax nmax];
amai = amai + amai0;

end

tia = toc; 
% disp(['Matriz do modelo feita em ', num2str(tia)])
 
Erro = 0;
catch
disp('Erro: Não foi possível criar a matriz do modelo')
Erro =1;
end
if Erro == 1
    continue
end

% Para gerar os limitantes na dimensão da matriz
[mt,nt] = size(PatK);
[mtd,ntd] = size(DemK);
dens = (numel([PatK; DemK]) - nnz([PatK; DemK]))*100/numel([PatK; DemK]);
amak = [amak; jt (mt+mtd) amai(1,1) nt amai(1,2) dens];

% informação da instância
datis0 = [jt k min(min(W)) max(max(W)) min(min(L)) max(max(L)) m min(w) max(w) min(l) max(l) min(d) max(d)];
datis = [datis; datis0];


% solucao relaxada
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);
tia = toc; 
disp(['Fez a matriz do modelo em ', num2str(tia)])
disp('Iniciando o CPLEX ...')

opt = cplexoptimset('cplex');
opt.timelimit = 3*60; 
% opt.display = 'on';

[x,fo,u,y] = cplexlp(ck,-DemK,-d,PatK,rhs,LB,UB,[],opt);
if isempty(fo)
    y
    disp('Erro: Cplex não resolveu o problema relaxado.')
    continue
end
I = find(ck);
% quantidade de objetos usados em cada tipo
x(I);
y;

% problema inteiro
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');
[xu,foi,u,yi] = cplexmilp(ck,-DemK,-d,PatK,rhs,[],[],[],LB,UB,type,[],opt);
if isempty(foi)
    yi
    continue
end
cys = [cys; [ num2str(jt),') ', yi.cplexstatusstring]];
% quantidade de objetos usados
obj = xu(I);
y;
% disp('Solucao inteira')
foI = foi;


% avalia algumas porcentagens
if ~isempty(xu)
    r = DemK*xu;
    tok = (L.*W)'*xu(I);
    pe = tok - ((r')*(l.*w));
    perda  = (pe*100)/tok;
    exce = (sum(r) - sum(d))*100;
    exc = exce/sum(d);
else    
    disp('ERROR')
    foi = 0;
    perda = 10000;
    exc = 10000;
end
pep = [jt perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];

n = length(I);

disp(' ' )
disp(' ' )
disp('Quantidade de objetos usados:')

for ku = 1:n
    if obj(ku,1) ~=0
      dok = [ku W(ku,1)  L(ku,1) obj(ku,1)];
      gei = [gei; dok];
    end
end

legdi = 'Objetos usados';
tit = 'Objeto & $W_{o}$ & $L_{o}$ & Quan. Uti. & Gcut';
LatexTab(gei,legdi,tit,[],1,[ ])
nam = num2str(jt);
csvwrite([csvdir,'/',nam,'WLusados.csv'],gei)

disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = 'Demanda e produção';
tit = 'Itens & $d_{i}$ & $Prod$';
ProtoP = [(1:mtd)' d r];
LatexTab(ProtoP,legdi,tit,[],1,[ ])
nam = num2str(jt);
csvwrite([csvdir,'/',nam,'DemPro.csv'],ProtoP)

% função objetivo, objetos e  tempo
fou0 = [jt fo foi tia y.time];
fou = [fou; fou0];
ortk = [jt obj(1,1) obj(2,1) obj(3,1) sum(obj)];
ort = [ort; ortk];
end

disp('========================================================' )
disp('RECORTAR AQUI')

% data
disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = 'Informações do problema';
tit = 'Prob. & o & W & L & m & w & $\ell$ & d';
LatexTab(datis,legdi,tit,[],1,[0 0 1 0 1 0 0 1 0 1 0 1 0])
csvwrite([csvdir,'/WLwld.csv'],datis)


% dimensao
disp(' ' )
disp(' ' )
disp('Dimensão da matriz:') 
legdi = 'Dimensão da matriz do modelo';
tit = 'Prob. & nr & lnr & nv & lnv & nz';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amak,legdi,tit,[],1,[])
csvwrite([csvdir,'/mnModelo.csv'],amak)


% objs
disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = 'Prob. & Obj. & Obj. & $o=1$ & $Go=1$ & $o=2$ & $Go=2$ & $o=3$ & $Go=3$ & ';
legdi = 'Dados da solução';
LatexTab(ort,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/Obj.csv'],ort)

% fou
disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = 'Prob. & Obj. & f.o.r & f.o.i & C & T. mod. & T. Cplex ';
legdi = 'Dados da solução';
LatexTab(fou,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/FOs.csv'],fou)

% perdas
disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = 'Análise da produção';
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/Perdas.csv'],peT)

cys

diary off