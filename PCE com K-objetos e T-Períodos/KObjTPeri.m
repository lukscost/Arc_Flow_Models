% Costa, L. L. S. (2022), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-estágios e com aparo.
% Programa para modelo de arco com múltiplos tipos de objetos e múltiplos
% períodos e que considera o estoque de itens

clc
clear all
close all
warning off 

diary('Output')
mkdir('CSV')

% car = 1
car = 2
if car == 1
disp('Sem custo de estoque e objeto vai ficando mais caro')
disp('$c_{z} = 1$')
disp('$c_{z} = (1 + 2*(t-1))$')
disp('$c_{y} = 0$')
else
disp('Valor inspirado na literatura')
pko = 0.001;
pki = 0.02;
disp(['$c_{z} = ',num2str(pko),'*L*W$'])
disp(['$c_{y} = ',num2str(pki),'*(w.*l)$'])
end

% tempo limite cplex
mik = 0.5;

opl

armaf = 1; 
peT =[];
funT = [];
dork =[];
amax = [];
cys = {};
objusa = []; 

for koi = 1:12

% dados
por = ['DadosTotal\KobjTpe',num2str(koi),'.mat'];
load(por);

disp(['Executando problema: ',num2str(koi),' ...'])

[w,I]=sort(w); l=l(I);
[md,nd] = size(d);

for i=1:nd
  d(:,i) = d(I,i);
end

% problema, W, L, quantidade de itens, periodos,
dri = [koi  min(W) max(W)  min(L) max(L) md min(w) max(w)  min(l) max(l)  min(min(d)) max(max(d))];
dork = [dork; dri];

disp('Iniciando montagem da matriz do modelo ...')

try
tic

[W,I]=sort(W); L=L(I); 

[mo,no] = size(W);

PatK0 = [];
DemK0 = [];
rhs0 = [];
ck0 = [];
amai0 = [0 0];
idek = [];
ieob = [];
% icy = zeros(1,mo);

for i=1:mo

id2es = [];

% preciso gerar as matrizes para os K objetos
[DE, U, Tt, gt, itTO, ft, I, o4, tew, cd, Dt, idz, id1es, idzs , id2es] = MatMode(L(i,1),l,W(i,1),w,d);

[mt,nt] = size(PatK0);
[mk,nk] = size(Tt);

if mt == 0 % PARA GERAR O BLOCO K=1
    PatK0 = Tt;
    DemK0 = DE;
    rhs0 = o4;

    % indice do vetor de custos
    inck = 1; 
else
    PatK0 = [PatK0 zeros(mt,nk); zeros(mk,nt) Tt];
    DemK0 = [DemK0  DE];
    rhs0 = [rhs0; o4];

    % indice do vetor de custos
    inck = [inck (nt+1)]; 
end

[mmax,nmax] = Ordem(W(i,1),L(i,1),l,w);
amaik = [mmax nmax];
amai0 = amai0 + amaik;

end

[mk,nk] = size(PatK0);
% [mrg,nrg] = size(rhs0)

% matriz para os itens guardados
Yes = sparse([sparse(zeros(mk,md)); sparse(eye(md))]);
Ver = sparse([sparse(PatK0); sparse(DemK0)]);
Pat = [Ver (-Yes)];

[mpa,npa] = size(Pat);
Dem = sparse([DemK0 zeros(md)]);

% indices do estoque de itens
inest = (1:nd)*npa;

PatK = [];
DemK = [];
rhs = [];
ck = [];
Ipro = [];
Iarma = [];
io = [];
IroK = [];

for i=1:nd

[mt,nt] = size(PatK);
    
if i == 1 % PARA GERAR O BLOCO i=1

    PatK = sparse(Pat);
    DemK = sparse([DemK Dem]);
    rhs = sparse([rhs0; d(:,i)]);
    c = zeros(npa,1);
    
    % CASO 1
if car == 1
    c(inck) = 1;    
    Iro = find(c);
    IroK = [IroK Iro']; 
else    
    % CASO 2
    ryu = (npa-md+1:npa);
    c(inck) = pko*(L.*W);
    Iro = find(c);
    IroK = [IroK Iro'];
    c(ryu) = pki*(w.*l);
end   
 
    ck = c;

else
    veh = sparse([sparse(mpa,(nt-md)) sparse(Yes)  sparse(Pat)]);

    PatK = sparse([PatK sparse(mt,(nk+md)); veh]);
    DemK = sparse([DemK Dem]);
    
    rhs = sparse([rhs; rhs0; d(:,i)]);
    
    [mt,nt] = size(PatK);
    
    c = zeros(npa,1);
    
    % CASO 1
if car == 1
    c(inck) = (1 + 2*(i-1));
    Iro = find(c);
    IroK = [IroK (Iro + (i-1)*npa)']; 
else    
     
    % CASO 2
    ryu = (npa-md+1:npa);
    c(inck) = pko*(L.*W);
    Iro = find(c);
    IroK = [IroK (Iro + (i-1)*npa)']; 
    c(ryu) = pki*(w.*l);
end          
    ck = [ck; c];
end
end
[mt,nt] = size(PatK);
IroK = [IroK nt];


% indice do objeto
Indp = inck;
for po = 2:nd
   Indp = [Indp (inck +(po-1)*npa)];
end

tia = toc; 
% disp(['Matriz do modelo feita em ', num2str(tia)])

mmax = amai0(1,1);
nmax = amai0(1,2);

Erro = 0;
catch
disp('Erro: Não foi possível criar a matriz do modelo')
Erro =1;
end

if Erro == 1
  amax0 =[koi inf nd*mmax inf nd*nmax inf];

    continue
end

dens = (numel(PatK) - nnz(PatK))*100/numel(PatK);
% (%) of the elements in the matrix are zeros.
amax0 =[koi mt  nd*mmax nt nd*nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

[mt,nt] = size(PatK);
 
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

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = mik*60; %8 mi
% opt.display = 'on';

[mck,nck] = size(ck);
[mpt,npt] = size(PatK);
[mrh,nrh] = size(rhs);

[x,fo,u,y] = cplexlp(ck,[],[],PatK,rhs,LB,UB,[],opt);
y;

% quantidade de objetos usados
I = find(ck);
x(I);
y;
foR = fo;

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');
% opt.display = 'on';

[xu,foi,u,yi] = cplexmilp(ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
cys = [cys; [ num2str(koi),') ', yi.cplexstatusstring]];
obj = xu(Indp);
yi;
foi;

WLa = [];
for uw=1:nd
WLa = [WLa (W.*L)'];
end

rh = length(IroK);
% avalia algumas porcentagens
if ~isempty(xu)
    r = DemK*xu;
    tok = WLa*xu(IroK(1:(rh-1)));
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
end

% indices PERIODOS E OBJETOS
pqi =[];
pqo =[];
for iu = 1:nd
    for iug = 1:mo
    pqi = [pqi iu];
    pqo = [pqo iug];
    end
end

% Produção nos períodos
produP = [];
for ij=1:nd
        inu = (1:npa) + (ij-1)*npa;
        produP  = [produP Dem*xu(inu)];
end
% produP

% produção de itens por objeto nos períodos
produO = [];
for ij=1:(rh-1)
    a = IroK(1,ij);
    b = (IroK(1,ij+1)-1);
    vrf = a:b;
    produO  = [produO DemK(:,vrf)*xu(vrf)];
end

% Quantidade de objetos estocados no periodo
esty =[];
for ieu = 1:nd
 iest = (inest(ieu)-md+1):(inest(ieu));
 esty = [esty xu(iest)];
end

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
ProtoP = [(1:md)' d(:,1) produP(:,1) d(:,2) produP(:,2)  d(:,3) produP(:,3) sum(d,2) r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(koi);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)

disp('                  ')

pep = [koi perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];

objpe = reshape(obj',[mo,nd]);
suobjk = sum(objpe,2);
suobjP = sum(objpe,1);

funp = [koi sum(suobjk)  suobjk' suobjP fo foi tia yi.time];
funT = [funT; funp];

disp(' ' )
disp('Uso dos objetos, t=1,2,3,...')
legdi = ['Uso dos objetos; Colunas: 1) Tipo do objeto; '...
    ' 2) Objetos usados no período 1'...
    ' 3) Objetos usados no período 2'...
    ' 4) Objetos usados no período 3'];
tit = 'Objeto & $t=1$ & $t=2$ & $t=3$';
Protok = [(1:mo)' objpe];
LatexTab(Protok,legdi,tit,[],1,[])
nam = num2str(koi);
csvwrite(['CSV/',nam,'DemPro.csv'],Protok)

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
LatexTab(dork,legdi,tit,[],1,[ ])
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
tit = ['Prob. & Obj. & $k=1$ & $k=2$ & $k=3$ & $t=1$ & $t=2$ & $t=3$'...
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
