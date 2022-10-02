% Costa, L. L. S. (2022), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-estágios e com aparo.
% Com múltiplos períodos e que considera estoque de itens e limite de
% estoque

clc
clear all
close all
warning off 

diary('Output')
mkdir('CSV')

% Casos
% car = 1
car = 1
if car == 1
disp('Sem custo de estoque e objeto vai ficando mais caro')
disp('$c_{z} = (1 + 2*(t-1))$')
% disp('$c_{z} = 1$')
disp('$c_{y} = 0$')
else
disp('Valor inspirado nos valores sugeridos na literatura')
disp('$c_{z} = 0.001*L*W$')
disp('$c_{y} = 0.02*(w.*l)$')
end

% limit and costs
limt = 500;
iotr = 50;
pko = 0.001;
pki = 0.02;

opl
armaf = 1; 

peT =[];
funT = [];
dork =[];
amax = [];
cys = {};

% time limit
mik = 0.5;

for koi = 1:12
por = ['Dados\Peri',num2str(koi),'.mat'];
load(por);

disp(['Executando problema: ',num2str(koi),' ...'])

[w,I]=sort(w); l=l(I);

[md,nd] = size(d);

for i=1:nd
  d(:,i) = d(I,i);
end

d;

% problema, W, L, quantidade de itens, periodos, 
dri = [koi W L md min(w) max(w)  min(l) max(l)  min(min(d)) max(max(d))];
dork = [dork; dri];

PatK = [];
DemK = [];
rhs = [];
ck = [];
Ipro = [];
Iarma = [];
%  pause
ikor = [];

disp('Iniciando montagem da matriz do modelo ...')

try
tic
% como o tamanho do objeto Ã© fixo sÃ³ preciso calcular uma vez
[DE, U, Tt, gt, itTO, ft, I, o4, tew, cd, Dt, idz, id1es, idzs , id2es] = MatMode(L,l,W,w,d);
[mk,nk] = size(Tt);

% matriz para os itens guardados
Yes = [zeros(mk,md); eye(md)];
Ver = [Tt; DE];
Pat = [Ver (-Yes)];
[mpa,npa] = size(Pat);
Dem = [DE zeros(md)];

for i=1:nd

[mt,nt] = size(PatK);
    
if i == 1 % PARA GERAR O BLOCO i=1

    PatK = Pat;
    DemK = [DemK Dem];
    
    Ip0 = (npa-md-cd+1):(npa-md);
    Ipro = [Ipro; Ip0];
    
    rhs = [o4; d(:,i)];
 
    c = zeros(npa,1);
    
    ikor = [ikor (npa-md+1:npa)'];
    % CASO 1
 if car == 1
    c(1) = 1;
else    
    % CASO 2
    ryu = (npa-md+1:npa);

    c(1) = pko*L*W;
    c(ryu) = pki*(w.*l);
end   
    ck = c;

else
    veh = [sparse(mpa,(nt-md)) Yes  Pat];

    PatK = [PatK sparse(mt,(nk+md)); veh];
    DemK = [DemK Dem];
    
    rhs = [rhs; o4; d(:,i)];
    
    Ip0 = (nt-md-cd+1):(nt-md);
    Ipro = [Ipro; Ip0];
    
    c = zeros(npa,1);
    
    [~,nt] = size(PatK);
    ikor = [ikor (nt-md+1:nt)'];
    
    % CASO 1
if car == 1
    c(1) = (1 + 2*(i-1));
else    
     
    % CASO 2
    ryu = (npa-md+1:npa);
     c(1) =  pko*L*W;
     c(ryu) = pki*(w.*l);
end
    ck = [ck; c];
end
end
[mt,nt] = size(PatK);
rag = zeros(nd,nt);
lid = limt*ones(nd,1);
inds = reshape(ikor,[md*nd,1]);

for caru = 1:nd
    rag(caru,ikor(:,caru)') = 1;
end

% indice do objeto
Indp = [];
for po = 1:nd
  Indp = [Indp (((po-1)*npa)+1)];
end

tia = toc; 
disp(['Matriz do modelo feita em ', num2str(tia)])

Erro = 0;
catch
disp('Erro: Não foi possível criar a matriz do modelo')
Erro =1;
end

if Erro == 1
    continue
end

[mt,nt] = size(PatK);
[mmax,nmax] = Ordem(W,L,l,w);

dens = (numel(PatK)-nnz(PatK))*100/numel(PatK);

amax0 =[koi mt  (nd*mmax) nt (nd*nmax) dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

% Permito armazenar item no fim ? 
% armaf = 0 NAO;
% armaf = 1 SIM
if armaf == 0
    PatK = PatK(:, (1:(nt-md)));
    DemK = DemK(:, (1:(nt-md)));
    ck = ck(1:(nt-md));
end

% Para gerar os limitantes
[mt,nt] = size(PatK);
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);
UB(inds) = iotr;

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = mik*60; %8 min

% com qtdade limite total de itens no esqtoque
[x,fo,u,y] = cplexlp(ck,rag,lid,PatK,rhs,LB,UB,[],opt);

% [x,fo,u,y] = cplexlp(ck,[],[],PatK,rhs,LB,UB,[],opt);
% cys = [cyz; y.cplexstatusstring]; 
if isempty(fo)
    y;
    disp('Erro: Cplex não resolveu o problema relaxado.')
    continue
end

% quantidade de objetos usados
I = find(ck);
x(I);
y;
foR = fo;

% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');
% opt.display = 'on';

% com qtdade limite total de itens no esqtoque
[xu,foi,u,yi] = cplexmilp(ck,rag,lid,PatK,rhs,[],[],[],LB,UB,type,[],opt);

% [xu,foi,u,yi] = cplexmilp(ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
cys = [cys; [ num2str(koi),') ', yi.cplexstatusstring]];

if isempty(foi)
    yi
    continue
end
% quantidade de objetos usados
obj = xu(Indp);
yi;
foi;

% avalia algumas porcentagens
if ~isempty(xu)
    r = DemK*xu;
    tok = (L*W)*xu(I);
    pe = sum(tok) - ((r')*(l.*w));
    perda  = (pe*100)/sum(tok);
    exce = (sum(r) - sum(sum(d)))*100;
    exc = exce/sum(sum(d));
    foi = foi;
    disp('Cplex encontrou uma solucao.')
else    
    disp('ERROR')
    foi = 0;
    perda = 10000;
    exc = 10000;
% end
% y
end


produ = zeros(md,nd);
for ij=1:nd
    produ(:,ij) = Dt*xu(Ipro(ij,:));
end

stoc = zeros(md,nd);
qtdst = zeros(1,nd);
for ik=1:nd
    xu(ikor(:,ik)');
    stoc(:,ik) = xu(ikor(:,ik)');
    qtdst(1,ik) = sum(stoc(:,ik));
end

disp(' ' )
disp('Itens, estoque, t=1,2,3,...')
legdi = ['Estoque; Colunas: 1) Itens; '...
    ' 2) Estoque do período 1'...
    ' 3) Estoque do período 2'...
    ' 4) Estoque do período 3'];
tit = 'Itens. &  $t=1$ &  $t=2$ &  $t=3$';
ProtES = [(0:md)' [qtdst; stoc]];
LatexTab(ProtES,legdi,tit,[],1,[])
nam = num2str(koi);
csvwrite(['CSV/',nam,'Estoq.csv'],ProtES)


% [W L];
disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = ['Demanda e produção; Colunas: 1) Itens; '...
    ' 2) Demanda do período 1'...
    ' 3) Produção do período 1'...
    ' 4) Demanda do período 2'...
    ' 5) Produção do período 2'...
    ' 6) Demadna do período 3'...
    ' 7) Produção do período 3'...
    ' 8) Demanda total'...
    ' 9) Produção total'];
tit = 'Itens. & $d_{1i}$ & $pro_{1i}$ & $d_{2i}$ & $pro_{2i}$ & $d_{3i}$ & $pro_{3i}$ & Dem. &  Prod';
ProtoP = [(1:md)' d(:,1) produ(:,1) d(:,2) produ(:,2)  d(:,3) produ(:,3) sum(d,2) r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(koi);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)

n = length(Indp);
disp('   ')
pep = [koi perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];

funp = [koi fo foi sum(obj) obj' tia yi.time];
funT = [funT; funp];
end

disp('========================================================' )
disp('RECORTAR AQUI')

disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = ['Informações do problema; Colunas: 1) Problema; '...
    ' 2) Valor da largura do objeto, W;' ...
    ' 3) Valor do comprimento do objeto, L;' ...
    ' 4) Quantidade de itens (m);'...
    ' 5) Intervalo da largura dos itens (w);'...
    ' 6) Intervalo do comprimento dos itens (l);'...
    ' 7) Intervalo da demanda de itens (d).'];
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 1 0 1 0 1 0])
csvwrite('CSV/WLwld.csv',dork)

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
LatexTab(amax,legdi,tit,[ ],1,[])
csvwrite('CSV/mnModelo.csv',amax)


disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = ['Prob. & Obj. & $t=1$ & $t=2$ & $t=3$'...
    '& f.o.r & f.o.i & T. mod. & T. Cplex '];
legdi = ['Dados da solução; Colunas: 1) Problema; '...
    ' 2) Total de objetos utilizados;' ...
    ' 3) Objetos usados no período $t=1$;' ...
    ' 4) Objetos usados no período $t=2$;' ...
    ' 5) Objetos usados no período $t=3$;' ...
    ' 6) Valor da função-objetivo do problema relaxado (f.o.r.);' ...
    ' 7) Valor da função-objetivo do problema inteiro (f.o.i);' ...
    ' 8) Tempo para construção da matriz do modelo;'...
    ' 9) Tempo para resolução do modelo.'];
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)


disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = ['Análise da produção; Colunas: 1) Problema; '...
    ' 2) Porcentagem de desperdício da área total dos objetos;' ...
    ' 3) Porcentagem do excesso de produção da quantidade total de itens produzidos;' ...
    ' 4) Quantidade de itens produzidos a mais.'];
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)

cys

diary off