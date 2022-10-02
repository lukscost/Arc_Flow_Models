% Costa, L. L. S. (2022), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-estágios e com aparo.
% Com múltiplos períodos e que considera estoque de itens e faixas (resíduos)

clc
clear all
close all
warning off 
opl

diary('Output')
mkdir('CSV')

pko = 0.001;
pki = 0.02;

% casos
car = 1
% car = 2
if car == 1
disp('Sem custo de estoque e objeto vai ficando mais caro')
% disp('$c_{z} = 1$')
disp('$c_{z} = (1 + 2*(t-1))$')
disp('$c_{q} = 0$')
disp('$c_{y} = 0$')
else
disp('Valor sugerido na literatura')
disp(['$c_{z} = ',num2str(pko),'*L*W$'])
disp(['$c_{q} = ',num2str(pko),'*L*ws$'])
disp(['$c_{y} = ',num2str(pki),'*(w.*l)$'])
end

armaf1 = 1; 
armaf2 = 1; 

dork = [];
cys = {};
amax = [];
funT = [];
peT = [ ];
badars = [];

% time limit
mik = 0.5;

for po = 1:12

% dados
por = ['Dados\Peri',num2str(po),'.mat'];
load(por);

disp(['Executando problema: ',num2str(po),' ...'])

[w,I]=sort(w); l=l(I);

[cdk,nd] = size(d);

for i=1:nd
  d(:,i) = d(I,i);
end

% problema, W, L, quantidade de itens, periodos, 
dri = [po W L cdk min(w) max(w)  min(l) max(l)  min(min(d)) max(max(d))];
dork = [dork; dri];

PatK = [];
DemK = [];
rhs = [];
ck = [];
Ipro = [];
Iarma = [];
io1 = [ ];
io2 = [ ];
r = [];
idz = [];
id1es = [];
idzs = [];
id2es = [];
idfa = [];
idites = [];
produ = [];
far =[];

disp('Iniciando montagem da matriz do modelo ...')

try
tic
% como o tamanho do objeto Ã© fixo sÃ³ preciso calcular uma vez
[DE, U, Tt, gt, itTO, ft, I, o4, tew, cd, Dt, Qtil, Ytil, idz, id1es, idzs , id2es] = MatMode(L,l,W,w,d);

% tamanho das matrizes
[lu,cu] = size(U);
[lde,cde] = size(DE);
[ldt,cdt] = size(Dt);
[mk,nk] = size(Tt);

idfa = (1:lu)+nk;
idites = (1:ldt)+nk+lu;

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

    % CASO 1
if car == 1
    c(1) = (1 + 2*(i-1));
    ryu2 = idfa; % faixas
    io2 = [io2 ryu2];
else
     
    % CASO 2
    ryu1 = idites; % itens
    ryu2 = idfa; % faixas
    c(1) = pko*L*W;
    c(ryu2) = pko*L*unique(w);
    c(ryu1) = pki*(w.*l);

    io1 = [io1 ryu1];
    io2 = [io2 ryu2];
end
            
    ck = c;
    %     indice do objeto
    Indp = idz;

else
    PatK = [PatK sparse(mt,npa); sparse(mpa,nt-nqtil-nytil) Qtil Ytil Pat];
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
%     c(1) = 1;
    ryu2 = idfa; % faixas
    io2 = [io2 ryu2 + (i-1)*npa];
else    
     
    % CASO 2
    ryu1 = idites; % itens
    ryu2 = idfa; % faixas

    c(1) = pko*L*W;
    c(ryu1) = pki*(w.*l);
    c(ryu2) = pko*L*unique(w);
    
    io1 = [io1 ryu1];
    io2 = [io2 ryu2 + (i-1)*npa];
end            
    ck = [ck; c];
    %     indice do objeto
    Indp = [Indp (((i-1)*npa)+1)];
end
end

[mt,nt] = size(PatK);

Indp;
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

% Permito armazenar itens no fim ? 
% armaf é definido no inicio
% armaf1 = 0 NAO;
% armaf1 = 1 SIM
if armaf1 == 0
    PatK = PatK(:, (1:(nt-lde)));
    DemK = DemK(:, (1:(nt-lde)));
    ck = ck(1:(nt-lde));
end

[mt,nt] = size(PatK);

% Permito armazenar faixas no fim ? 
% armaf = 0 NAO;
% armaf = 1 SIM
if armaf2 == 0
    PatK = PatK(:, (1:(nt-lu)));
    DemK = DemK(:, (1:(nt-lu)));
    ck = ck(1:(nt-lu));
end

[mmax,nmax] = Ordem(W,L,l,w);
% (%) of the elements in the matrix are zeros.
dens = (numel(PatK)-nnz(PatK))*100/numel(PatK);
amax0 =[po mt  nd*mmax nt nd*nmax dens];
amax = [amax; amax0];

disp('Iniciando o CPLEX ...')

% Para gerar os limitantes
LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = mik*60; %8 min

[x,fo,u,y] = cplexlp(ck,[],[],PatK,rhs,LB,UB,[],opt);
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
fo;

% disp('Solucao inteira')

% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');
% opt.display = 'on';

[xu,foi,u,yi] = cplexmilp(ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
cys = [cys; [ num2str(po),') ', yi.cplexstatusstring]];
if isempty(foi)
    yi;
    disp('Erro: Cplex não resolveu o problema inteiro.')
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
% end
% y
end

% produção de itens nos períodos
for i=1:nd
    pef = id2es + (i-1)*npa;
    produ(:,i) = Dt*xu(pef,1);
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
ProtoP = [(1:cdk)' d(:,1) produ(:,1) d(:,2) produ(:,2)  d(:,3) produ(:,3) sum(d,2) sum(r,2)];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(po);
csvwrite(['CSV/',nam,'DemPro.csv'],ProtoP)
disp(' ' )
disp(' ' )

Estoqs = [(1:cdk)' (produ(:,1)-d(:,1)) (produ(:,1)-d(:,1)+produ(:,2)-d(:,2))  (produ(:,1)-d(:,1)+produ(:,2)-d(:,2)+produ(:,3)-d(:,3))];
legdi = 'Armazenamento de itens nos períodos';
tit = '# & 1 & 2 & 3';
LatexTab(Estoqs,legdi,tit,[],1,[])
csvwrite(['CSV/estoqus',num2str(po),'.csv'],Estoqs)

funp = [po sum(obj) obj' fo foi tia yi.time];
funT = [funT; funp];

pep = [po perda exc  (sum(r) - sum(sum(d)))];
peT = [peT; pep];

n = length(Indp);

disp(' ' )
disp(' ' )
disp('Quantidade de faixas guardadas no último periodo:')
vers = [(1:length(unique(w)))' reshape(xu(io2),[length(unique(w)),3])];

legdi = 'Armazenamento de faixas nos períodos';
tit = '# & 1 & 2 & 3';
LatexTab(vers,legdi,tit,[],1,[])
csvwrite(['CSV/vers',num2str(po),'.csv'],vers)
disp(' ' )

badars = [badars; sum(vers(:,2:4).*unique(w))];

[mt,nt] = size(PatK);
far = xu((nt-nytil-nqtil+1):(nt-nytil));
wut = unique(w);
if armaf2 == 0
    disp('Não foi permitido guardar faixas no último periodo.')
else
    if isempty(find(far,1)) && (armaf2 == 1)
        disp('Não foi armazenada nenhuma faixa no último periodo.')
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

badars

cys
diary off