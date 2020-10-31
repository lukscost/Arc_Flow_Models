% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios,
% Modelo de fluxo em arcos tradicional

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

% "for" para resolver múltiplos períodos individualmente    
% for kio = 1:3
kio = 1;
    
% dados
por = ['Gcut\gcut',num2str(po),'d.mat'];
load(por)
[w,I]=sort(w); l = l(I); dr = d(I,:); 
d = dr(:,kio);
% [w l d];

% quantidade de itens, sem rotação e com 
mr = length(w);
% indice dos itens
ri = (1:mr)';

if kio == 1
% problema, W, L, itens, demanda 
dri = [po W L mr min(w) max(w)  min(l) max(l)  min(d) max(d)];
dork = [dork; dri];

% Matriz do modelo
try
tic
disp('Iniciando montagem da matriz do modelo ...')
[DE, U, Tt, gt, itTO, ft, o4] = MatMode(L,l,W,w,d,ri,mr);

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
end

% Para gerar os limitantes na dimensão da matriz
[mt,nt] = size(Tt);
[m1,n1] = size([Tt; DE]);
[mmax,nmax] = Ordem(W,L,l,w,w);
dens = nnz([Tt; DE]==0)*100/numel([Tt; DE]);
if kio == 1
amax0 =[po m1 mmax n1 nmax dens];
amax = [amax; amax0];
end


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


fok(kio,1)  = fo;
foik(kio,1) = foi;
tiak(kio,1) =  tia;
yit(kio,1) = yi.time ;

% if kio == 3
funp = [po sum(foik) fok foik' mean(tia) mean(yit)];
funT = [funT; funp];
% end

% avalia algumas porcentagens
if ~isempty(xu)
    r = DE*xu;
    rou(kio,1) = (r')*(l.*w);
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
end

% informações das perdas
pep = [po perda exc  (sum(r) - sum(d))];
peT = [peT; pep];

% [W L]
% [w l d r]
disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = 'Demanda e produção';
tit = 'Itens & demanda-i & produção-i,';
ProtoP = [(1:mr)' d  r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(po);
% end

% avalia algumas porcentagens
aryt = sum(peT(:,4));
excer = (aryt - sum(sum(dr)))*100;
excr = exce/sum(d);
tokl = L*W*sum(foik);
per = tokl - sum(rou);
perk  = (per*100)/tokl;
perg(po,:) = [po perk excr aryt];
peT=[];

end


disp('========================================================' )
disp('RECORTAR AQUI')

disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = 'Informações do problema';
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 1 0 1 0 1 0])
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
LatexTab(perg,legdi,tit,[ ],1,[])
csvwrite('CSV/Perdas.csv',peT)

disp(cys)

diary off