% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios,
% com aparo e múltiplos períodos e que permite estoque de itens

clc
clear all
close all
warning off 

% script para salvar o caminho das funções do CPLEX
opl

diary('Testes')
mkdir('CSV')

% custos da função objetivo
pko = 0.001;
pki = 0.02;

% CASOS
car = 1
% car = 2
if car == 1
disp('Sem custo de estoque e objeto vai ficando mais caro')
disp('$c_{z} = (1 + 2*(t-1))$')
disp('$c_{y} = 0$')
else
disp('Valor inspirado nos valores sugeridos na literatura')
disp(['$c_{z} = ',num2str(pko),'*L*W$'])
disp(['$c_{y} = ',num2str(pki),'*(w.*l)$'])
end

% armazena itens ou faixas no final ? Sim = 1
armaf = 1; 

% matrizes dos dados
peT =[]; funT = []; dork =[]; 
amax = []; cys = {};


for koi = 1:1
    
% load dos dados
por = ['Dados\Peri',num2str(koi),'.mat'];
load(por);

disp(['Executando problema: ',num2str(koi),' ...'])

% organiza os dados
[w,I]=sort(w); l=l(I);
[md,nd] = size(d);
for i=1:nd
  d(:,i) = d(I,i);
end


% problema, W, L, quantidade de itens, periodos, 
dri = [koi W L md min(w) max(w)  min(l) max(l)  min(min(d)) max(max(d))];
dork = [dork; dri];

% inicia as matrizes do modelo
PatK = []; DemK = []; rhs = [];
ck = []; Ipro = []; Iarma = []; io = [];

disp('Iniciando montagem da matriz do modelo ...')


tic
try
% função responsável pela criação das matrizes de cada tipo de objeto
% como o tamanho do objeto são fixo só precisamos calcular uma vez
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
    
    % CASO 1
 if car == 1
    c(1) = 1;
else    
    % CASO 2
    ryu = (npa-md+1:npa);
    c(1) = pko*L*W;
    c(ryu) = pki*(w.*l);
    iesto = npa-md+1:npa;
    io = [io iesto];
end   
    ck = c;

else
    veh = [zeros(mpa,(nt-md)) Yes  Pat];
    PatK = [PatK zeros(mt,(nk+md)); veh];
    DemK = [DemK Dem];
    rhs = [rhs; o4; d(:,i)];
    [mt,nt] = size(PatK);
    Ip0 = (nt-md-cd+1):(nt-md);
    Ipro = [Ipro; Ip0];
    c = zeros(npa,1);
    
    % CASO 1
if car == 1
    c(1) = (1 + 2*(i-1));
else    
     
    % CASO 2
    ryu = (npa-md+1:npa);
    c(1) =  pko*L*W;
    c(ryu) = pki*(w.*l);
    iesto = nt-md+1:nt;
    io = [io iesto];
end
    ck = [ck; c];

end
end

% indice do objeto
Indp = [];
for po = 1:nd
  Indp = [Indp (((po-1)*npa)+1)];
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
[mmax,nmax] = Ordem(W,L,l,w);
dens = (numel(PatK)-nnz(PatK))*100/numel(PatK);
amax0 =[koi mt  (nd*mmax) nt (nd*nmax) dens];
amax = [amax; amax0];


% Geralmente a retirada dessas partes deixa o modelo infactível.
% Permito armazenar itens e faixas no fim ? 
% armaf = 0 NAO;
% armaf = 1 SIM
if armaf == 0
    PatK = PatK(:, (1:(nt-md)));
    DemK = DemK(:, (1:(nt-md)));
    ck = ck(1:(nt-md));
end

disp('Iniciando o CPLEX ...')

% Para gerar os limitantes
[mt,nt] = size(PatK);
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = 4*60; %4 min
[x,fo,u,y] = cplexlp(ck,[],[],PatK,rhs,LB,UB,[],opt);
if isempty(fo)
    y
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
[xu,foi,u,yi] = cplexmilp(ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
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
end


fari = xu((nt-md+1):nt);
if isempty(find(fari,1)) && (armaf == 1)
     disp('Não foi armazenado nenhum item no último periodo.')
else
     disp('Itens armazenados:')
     LatexTab(fari','Itens armazenados',num2str(1:md),[],1,[])
end


produ = [];
for ij=1:nd
    produ(:,ij) = Dt*xu(Ipro(ij,:));
end

% [W L];
disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = 'Demanda e produção';
tit = 'Itens. & $d_{1i}$ & $pro_{1i}$ & $d_{2i}$ & $pro_{2i}$ & $d_{3i}$ & $pro_{3i}$ & Dem. &  Prod';
ProtoP = [(1:md)' d(:,1) produ(:,1) d(:,2) produ(:,2)  d(:,3) produ(:,3) sum(d,2) r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(koi);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)

n = length(Indp);

disp('                  ')

% armazena os dados dos resultados
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
% Co = num2str((1:5)');
% gol = native2unicode(double(unicode2native(')'))*(ones(5,1)));
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
% Co = num2str((1:5)');
% gol = native2unicode(double(unicode2native(')'))*(ones(5,1)));
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
% Co = num2str((1:4)');
% gol = native2unicode(double(unicode2native(')'))*(ones(4,1)));
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)

cys

diary off
save('All.mat')