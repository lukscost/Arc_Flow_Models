% Costa, L. L. S. (2021), Extens�es do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matem�tica Estat�stica e Computa��o Cient�fica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-est�gios,
% com aparo e m�ltiplos per�odos e que permite estoque de itens e faixas

clc
clear all
close all
warning off 

% script para salvar o caminho das fun��es do CPLEX
opl

diary('Testes')
mkdir('CSV')

% custos da fun��o objetivo
pko = 0.001;
pki = 0.02;

% CASOS
car = 1
% car = 2
if car == 1
disp('Sem custo de estoque e objeto vai ficando mais caro')
disp('$c_{z} = (1 + 2*(t-1))$')
disp('$c_{q} = 0$')
disp('$c_{y} = 0$')
else
disp('Valor sugerido na literatura')
disp(['$c_{z} = ',num2str(pko),'*L*W$'])
disp(['$c_{q} = ',num2str(pko),'*L*ws$'])
disp(['$c_{y} = ',num2str(pki),'*(w.*l)$'])
end


% armazena itens ou faixas no final ? Sim = 1
armaf1 = 1; 

% matrizes dos dados
dork = []; cys = {};
amax = []; funT = [];
peT = [ ];


for po = 1:1
    
% load dos dados
por = ['Dados\Peri',num2str(po),'.mat'];
load(por);

disp(['Executando problema: ',num2str(po),' ...'])

% organiza os dados
[w,I]=sort(w); l=l(I);
[cdk,nd] = size(d);
for i=1:nd
  d(:,i) = d(I,i);
end
d;


% problema, W, L, quantidade de itens, periodos, 
dri = [po W L cdk min(w) max(w)  min(l) max(l)  min(min(d)) max(max(d))];
dork = [dork; dri];


% inicia as matrizes do modelo
PatK = []; DemK = []; rhs = []; ck = []; Ipro = []; Iarma = [];
io1 = [ ]; io2 = [ ]; r = []; idz = []; id1es = [];
idzs = []; id2es = []; idfa = []; idites = []; produ = []; far =[];

disp('Iniciando montagem da matriz do modelo ...')

try
tic
% fun��o respons�vel pela cria��o das matrizes de cada tipo de objeto
% como o tamanho do objeto s�o fixo s� precisamos calcular uma vez
[DE, U, Tt, gt, itTO, ft, I, o4, tew, cd, Dt, Qtil, Ytil, idz, id1es, idzs , id2es] = MatMode(L,l,W,w,d);


% tamanho das matrizes
[lu,cu] = size(U);
[lde,cde] = size(DE);
[ldt,cdt] = size(Dt);
[mk,nk] = size(Tt);
idfa = (1:lu)+nk;
idites = (1:ldt)+nk+lu;
% blocos
Pat = [[Tt; DE] -Qtil -Ytil];
[mpa,npa] = size(Pat);
Dem = [DE zeros(lde,lu+lde)];
[ldem,cdem] = size(Dem);
[mqtil, nqtil] = size(Qtil);
[mytil, nytil] = size(Ytil);


for i=1:nd

[mt,nt] = size(PatK);
    
if i == 1 % PARA GERAR O BLOCO i=1

    PatK = Pat;
    DemK = [DemK Dem];

    rhs = [rhs; o4; d(:,i)];

    Ip0 = (npa-cdt-lu+1):(npa-cdt-lu);
    Ipro = [Ipro; Ip0];
    
    % vetor de custos do PI
    c = zeros(npa,1);

if car == 1
    c(1) = (1 + 2*(i-1));
else
     
    % CASO 2
    ryu1 = idites; % itens
    ryu2 = idfa; % faixas
    c(1) = pko*L*W;
    c(ryu1) = pko*(w.*l);
    c(ryu2) = pki*L*unique(w);
    io1 = [io1 ryu1];
    io2 = [io2 ryu2];
end
            
    ck = c;
    %     indice do objeto
    Indp = idz;

else
    PatK = [PatK zeros(mt,npa);zeros(mpa,nt-nqtil-nytil) Qtil Ytil Pat];
    DemK = [DemK Dem];
    rhs = [rhs; o4; d(:,i)];
    [mt,nt] = size(PatK);
    
    Ip0 = (nt-lu-cdt+1):(nt-cdt-lu);
    Ipro = [Ipro; Ip0];
    
    % vetor de custos do PI
    c = zeros(npa,1);

    % CASO 1
if car == 1
     c(1) = (1 + 2*(i-1));
else    
    % CASO 2
    ryu1 = idites; % itens
    ryu2 = idfa; % faixas
    c(1) = pko*L*W;
    c(ryu1) = pko*(w.*l);
    c(ryu2) = pki*L*unique(w);
    io1 = [io1 ryu1];
    io2 = [io2 ryu2];
end            
    ck = [ck; c];
    %     indice do objeto
    Indp = [Indp (((i-1)*npa)+1)];
end
end

[mt,nt] = size(PatK);
Indp;
tia = toc; 
% disp(['Matriz do modelo feita em ', num2str(tia)])
Erro = 0;
catch
disp('Erro: N�o foi poss�vel criar a matriz do modelo')
Erro =1;
end
if Erro == 1
    continue
end

% Geralmente a retirada dessas partes deixa o modelo infact�vel.
% Permito armazenar itens e faixas no fim ? 
% armaf = 0 NAO;
% armaf = 1 SIM
[mt,nt] = size(PatK);
if armaf1 == 0
    PatK = PatK(:, (1:(nt-lde - lu)));
    DemK = DemK(:, (1:(nt-lde - lu)));
    ck = ck(1:(nt-lde- lu));
end

% Para gerar os limitantes na dimens�o da matriz
[mmax,nmax] = Ordem(W,L,l,w);
dens = (numel(PatK)-nnz(PatK))*100/numel(PatK);
amax0 =[po mt  nd*mmax nt nd*nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

% disp('Solucao relaxada')
% Para gerar os limitantes
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);
opt = cplexoptimset('cplex');
opt.timelimit = 1*60; %4 min

% CPLEX
[x,fo,u,y] = cplexlp(ck,[],[],PatK,rhs,LB,UB,[],opt);
if isempty(fo)
    y;
    disp('Erro: Cplex n�o resolveu o problema relaxado.')
    continue
end
% quantidade de objetos usados
I = find(ck);
x(I);
y;
fo;

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');

[xu,foi,u,yi] = cplexmilp(ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
cys = [cys; yi.cplexstatusstring]; 
if isempty(foi)
    yi;
    disp('Erro: Cplex n�o resolveu o problema inteiro.')
    continue
end
% quantidade de objetos usados
yi;
foi;
obj = xu(Indp);


% avalia algumas porcentagens
if ~isempty(xu)
    r = DemK*xu;
    tok = (L*W)*sum(xu(I));
    pe = sum(tok) - ((sum(r,2)')*(l.*w));
    perda  = (pe*100)/sum(tok);
    exce = (sum(sum(r)) - sum(sum(d)))*100;
    exc = exce/sum(sum(d));
    fou = foi;
else    
    disp('ERROR')
    fou  = 0;
    perda = 10000;
    exc  = 10000;
end

% produ��o de itens nos per�odos
for i=1:nd
    pef = id2es + (i-1)*npa;
    produ(:,i) = Dt*xu(pef,1);
end


% [W L];
disp(' ' )
disp('Itens, demanda-i, produ��o-i, i=1,2,3')
legdi = 'Demanda e produ��o';
tit = 'Itens. & $d_{1i}$ & $pro_{1i}$ & $d_{2i}$ & $pro_{2i}$ & $d_{3i}$ & $pro_{3i}$ & Dem. &  Prod';
ProtoP = [(1:cdk)' d(:,1) produ(:,1) d(:,2) produ(:,2)  d(:,3) produ(:,3) sum(d,2) sum(r,2)];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(po);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)



% armazena os dados dos resultados
funp = [po sum(obj) obj' fo foi tia yi.time];
funT = [funT; funp];
pep = [po perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];
n = length(Indp);

disp(' ' )
disp(' ' )
disp('Quantidade guardada no �ltimo periodo:')
[mt,nt] = size(PatK);
far = xu((nt-nytil-nqtil+1):(nt-nytil));
fari = xu((nt-nytil+1):nt);

wut = unique(w);
if armaf1 == 0
    disp('N�o foi permitido guardar faixas no �ltimo periodo.')
else
    if isempty(find(fari,1)) && (armaf1 == 1)
            disp('N�o foi armazenado nenhum item no �ltimo periodo.')
    else
        disp('Itens armazenados:')
        LatexTab(fari','Itens armazenados',num2str(1:nytil),[],1,[])
    end
     
    
    if isempty(find(far,1)) && (armaf1 == 1)
        disp('N�o foi armazenada nenhuma faixa no �ltimo periodo.')
    else
        for k = 1:length(wut)
            if far(k,1) ~= 0
               Irow = find(w==wut(k));
               if length(Irow) >= 2
     arew = ['Foram armazenadas ', num2str(far(k,1)), ' faixas de largura '...
    ,num2str(wut(k)),'x',num2str(L),', largura do item: '];
for hou = 1:length(Irow)
    if hou ~=length(Irow)
    arew = [arew, num2str(Irow(hou)), ' e '];
    else
    arew = [arew, num2str(Irow(hou)), '.'];
    end
end           
disp(arew)
               else
disp(['Foram armazenadas ', num2str(far(k,1)), ' faixas de largura '...
    ,num2str(wut(k)),'x',num2str(L),', largura do item: '...
    ,num2str(find(w==wut(k)))])
               end
            end
        end
    end
end

end

disp('========================================================' )
disp('RECORTAR AQUI')

disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = 'Informa��es do problema';
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 1 0 1 0 1 0])
csvwrite('CSV/WLwld.csv',dork)


disp(' ' )
disp(' ' )
disp('Dimens�o da matriz:') 
legdi = 'Dimens�o da matriz do modelo';
tit = 'Prob. & nr & lnr & nv & lnv & nz';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amax,legdi,tit,[ ],1,[])
csvwrite('CSV/mnModelo.csv',amax)


disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = ['Prob. & Obj. & $t=1$ & $t=2$ & $t=3$'...
    '& f.o.r & f.o.i & T. mod. & T. Cplex '];
legdi = 'Dados da solu��o';
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSV/FOs.csv',funT)


disp(' ' )
disp(' ' )
disp('desperdicio da �rea total \% , exec produ \%, exec em itens') 
legdi = 'An�lise da produ��o';
tit = 'Prob. & $(\%$ - desperd�cio$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)


cys
diary off