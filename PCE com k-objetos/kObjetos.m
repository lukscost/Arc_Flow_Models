% Costa, L. L. S. (2022), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-estágios e com aparo.
% Com múltiplos tipos de objetos

clc
clear all
close all
warning off 
format long

diary('Output')
csvdir = 'CSVGcut';
mkdir(csvdir)

disp('Ax >= d')

% limit time to find solution (min)
mik = 0.5;

% display log of CPLEX ? 1 = Yes
dispsh = 0;

opl

fou = [];
datis = [];
amak = [];
dork = [];
cys = {};
peT = [];
ort = [];


for jt = 1:12

gei = [];

disp(['Executando problema: ',num2str(jt),' ...'])
runjt = ['DadosGcut/gcut',num2str(jt),'d.m'];
run(runjt)

k = length(W);
m = length(w);

disp(['Dimensão: ',num2str(k),'x',num2str(m)])

%  organização dos dados
[w,I]=sort(w); l=l(I); d = d(I);
[mo,no] = size(W);

PatK = [];
DemK = [];
rhs = [];
ck = [];
amai = [0 0];

disp('Iniciando montagem da matriz do modelo ...')

try
tic;

for i=1:mo
     
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
    PatK = [PatK sparse(mt,nk); sparse(mk,nt) Tt];
    DemK = [DemK  DE];
    rhs = [rhs; o4];
    
    % vetor de custos do PI
    c = zeros(nk,1);
    c(1) = L(i,1)*W(i,1);
    ck = [ck; c];
end

[mmax,nmax] = Ordem(W(i,1),L(i,1),l,w);
amai0 = [mmax nmax];
amai = amai + amai0;

end
 
Erro = 0;
catch
disp('Erro: Não foi possível criar a matriz do modelo')
Erro =1;
end

if Erro == 1
    continue
end

% Para gerar os limitantes
[mt,nt] = size(PatK);
[mtd,ntd] = size(DemK);

datis0 = [jt k min(min(W)) max(max(W)) min(min(L)) max(max(L)) m min(w) max(w) min(l) max(l) min(d) max(d)];
datis = [datis; datis0];

dens = (numel([PatK; DemK]) - nnz([PatK; DemK]))*100/numel([PatK; DemK]);
amak = [amak; jt (mt+mtd) amai(1,1) nt amai(1,2) dens];

LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);

tia = toc; 
disp(['Fez a matriz do modelo em ', num2str(tia)])
disp('Iniciando o CPLEX ...')

% tipo da solucao relaxada
opt = cplexoptimset('cplex');
opt.timelimit = mik*60; %8 min

[x,fo,u,y] = cplexlp(ck,-DemK,-d,PatK,rhs,LB,UB,[],opt);

if isempty(fo)
    y
    disp('Erro: Cplex não resolveu o problema relaxado.')
    continue
end

% quantidade de objetos usados
I = find(ck);
x(I);
y;
% disp('Solucao relaxada')
% fou0 = fo;

% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');
 
if dispsh == 1
 opt.display = 'on';
else
 opt.display = 'off';
end

[xu,foi,u,yi] = cplexmilp(ck,-DemK,-d,PatK,rhs,[],[],[],LB,UB,type,[],opt);

% quantidade de objetos usados
cys = [cys; [ num2str(jt),') ', yi.cplexstatusstring]];

if isempty(foi)
    yi
    continue
end
obj = xu(I);
y;
% disp('Solucao inteira')
% foI = foi;

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
% end
% y
end

pep = [jt perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];

n = length(I);

disp(' ' )
disp(' ' )
disp('Quantidade de objetos usados:')

for ku = 1:n
    if obj(ku,1) ~=0
%         disp(['Objeto ',num2str(k), ' (',num2str(L(k,1)),'x',num2str(W(k,1)),...
%             '): foi usado ', num2str(obj(k,1)),' vezes.'])
      dok = [ku W(ku,1)  L(ku,1) obj(ku,1)];
      gei = [gei; dok];
    end
end

% gei
legdi = ['Objetos usados; Colunas: 1) Número do objeto; '...
    ' 2) Largura do objeto $W_{o}$; '...
    ' 2) Comprimento do objeto $L_{o}$; '...
    ' 3) Quantidade utilizada; 4) Quantidade utilizada GCUT.'];
tit = 'Objeto & $W_{o}$ & $L_{o}$ & Quan. Uti. & Gcut';
LatexTab(gei,legdi,tit,[],1,[ ])
nam = num2str(jt);
csvwrite([csvdir,'/',nam,'WLusados.csv'],gei)

disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = ['Demanda e produção; Colunas: 1) Itens; '...
    ' 2) Demanda do período;'...
    ' 3) Produção do período.'];
tit = 'Itens & $d_{i}$ & $Prod$';
ProtoP = [(1:mtd)' d r];
LatexTab(ProtoP,legdi,tit,[],1,[ ])
nam = num2str(jt);
csvwrite([csvdir,'/',nam,'DemPro.csv'],ProtoP)

fou0 = [jt fo foi tia y.time];
fou = [fou; fou0];

ortk = [jt obj(1,1) obj(2,1) obj(3,1) sum(obj)];
ort = [ort; ortk];

end

disp('========================================================' )
disp('RECORTAR AQUI')

% datis
disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = ['Informações do problema; Colunas: 1) Problema; '...
    ' 2) Quantidade de objetos (o);'...
    ' 3) Largura do objeto, W;' ...
    ' 4) Comprimento do objeto, L;' ...
    ' 5) Quantidade de itens (m);'...
    ' 6) Intervalo da largura dos itens (w);'...
    ' 7) Intervalo do comprimento dos itens (l);'...
    ' 8) Intervalo da demanda de itens (d).'];
tit = 'Prob. & o & W & L & m & w & $\ell$ & d';
LatexTab(datis,legdi,tit,[],1,[0 0 1 0 1 0 0 1 0 1 0 1 0])
csvwrite([csvdir,'/WLwld.csv'],datis)


% dimensao
disp(' ' )
disp(' ' )
disp('Dimensão da matriz:') 
legdi = ['Dimensão da matriz do modelo; Colunas: 1) Problema; '...
    ' 2) Número de restrições (nr);' ...
    ' 3) Limitante no número de restrições (lnr);' ...
    ' 4) Número de variáveis; (nv);' ...
    ' 5) Limitante no número de variáveis (lnv);'...
    ' 6) $\%$ de elementos nulos (nz).'];
tit = 'Prob. & nr & lnr & nv & lnv & nz';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amak,legdi,tit,[],1,[])
csvwrite([csvdir,'/mnModelo.csv'],amak)


% fou
disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = ['Prob. & Obj. & Obj. & $o=1$ & $Go=1$ & $o=2$ & $Go=2$ & $o=3$ & $Go=3$ & '];
legdi = ['Dados da solução; Colunas: 1) Problema; '...
    ' 2) Total de objetos utilizados;' ...
    ' 3) Total de objetos utilizados (Gcut);' ...
    ' 4) Objetos utilizados do tipo $o=1$;' ...
    ' 5) Objetos utilizados do tipo $o=1$ (Gcut);' ...
    ' 6) Objetos utilizados do tipo $o=2$;' ...
    ' 7) Objetos utilizados do tipo $o=2$ (Gcut);' ...
    ' 8) Objetos utilizados do tipo $o=3$;' ...
    ' 9) Objetos utilizados do tipo $o=3$ (Gcut).'];
LatexTab(ort,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/FOs.csv'],ort)

% fou
disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = ['Prob. & Obj. & f.o.r & f.o.i & C & T. mod. & T. Cplex '];
legdi = ['Dados da solução; Colunas: 1) Problema; '...
    ' 2) Total de objetos utilizados;' ...
    ' 3) Objetos utilizados do tipo $o=1$;' ...
    ' 4) Objetos utilizados do tipo $o=2$;' ...
    ' 5) Objetos utilizados do tipo $o=3$;' ...
    ' 6) Valor da função-objetivo do problema relaxado (f.o.r.);' ...
    ' 7) Valor da função-objetivo do problema inteiro (f.o.i);' ...
    ' 8) Valor da função-objetivo (Cintra);' ...
    ' 9) Tempo para construção da matriz do modelo;'...
    ' 10) Tempo para resolução do modelo.'];
LatexTab(fou,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/FOs.csv'],fou)

disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = ['Análise da produção; Colunas: 1) Problema; '...
    ' 2) Porcentagem de desperdício da área total dos objetos;' ...
    ' 3) Porcentagem do excesso de produção da quantidade total de itens produzidos;' ...
    ' 4) Quantidade de itens produzidos a mais.'];
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/Perdas.csv'],peT)

cys

diary off