% Costa, L. L. S. (2021), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do PCE bidimensional, guilhotinado, 2-estágios.
% Decomposição do modelo de fluxo em arcos e que permite rotação de 
% 90º nos itens.

clc
clear all
close all
warning off 
diary('TestesCOM')

mkdir('CSVcom')
disp('Decomposição COM rotação')

% script para salvar o caminho das funções do CPLEX
opl

% matrizes dos dados
dork = []; amax = []; funT = [];
peT =  []; cys = {};


for po = 1:1

% dados
por = ['Gcut\gcut',num2str(po),'d.mat'];
load(por);

[w,I]=sort(w); l = l(I); d = d(I); 
% [w l];

disp(['Executando problema: ',num2str(po),' ...'])

% quantidade de itens, sem rotação e com 
mwf = length(w);

% organiza os dados
% para nao duplicar itens quadrados
Id = find((w-l)==0);
rota = (1:mwf)';
rota(Id) = [];
% para rotacionar os itens que nao sao quadrados
ri = [(1:mwf)'; rota];
wr = [w; l(rota)];
lr = [l; w(rota)];
dr = [d; d(rota)];
mr = length(wr);
mu = length(unique(wr));
[wr,I]=sort(wr); lr = lr(I);  ri = ri(I);  dr = dr(I); 


% problema, W, L, itens, demanda 
dri = [po W L mwf min(w) max(w)  min(l) max(l)  min(d) max(d)];
dork = [dork; dri];

% Matriz do modelo
try
tic
disp('Iniciando montagem da matriz do modelo ...')
[A1eq, A1de, b1eq, A2eq, A2de, b2eq] = MatMode(L,lr,W,wr,ri,mwf);
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


disp('Iniciando o CPLEX ...')

% CORTE DOS ITENS
[m2t,n2t] = size(A2de);
[m2tq,n2tq] = size(A2eq);

% vetor de custos do PI
c = zeros(n2t,1);
c(1:mu,1) = unique(wr);

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(n2t,1);
type = char(tipe');

% Para gerar os limitantes
LB = sparse(zeros(n2t,1));
UB = inf*ones(n2t,1);

% disp('Solucao relaxada')
opt = cplexoptimset('cplex');
opt.timelimit = 8*60; %8 min
[xu1,fo1,u1,yi1] = cplexmilp(c,-A2de,-d,A2eq,b2eq,[],[],[],LB,UB,type,[],opt);
yu = xu1(1:mu);

% CORTE DAS FAIXAS
% tamanho das matrizes
[m1t,n1t] = size(A1de);
[m1tq,n1tq] = size(A1eq);

% Para gerar os limitantes na dimensão da matriz
% [mmax,nmax] = Ordem(W,L,lr,wr,d)
dens1 = (numel([A1de; A1eq])-nnz([A1de; A1eq]))*100/numel([A1de; A1eq]);
dens2 = (numel([A2de; A2eq])-nnz([A2de; A2eq]))*100/numel([A2de; A2eq]);
amax0 =[po (m2t+m2tq)  n2t dens2 (m1t+m1tq)  n1t dens1 (m1t+m1tq+m2t+m2tq) (n2t+n1t)];
amax = [amax; amax0];

% vetor de custos do PI
c = zeros(n1t,1);
c(1) = 1;

% disp('Solucao inteira')
% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(n1t,1);
type = char(tipe');

% Para gerar os limitantes
LB = sparse(zeros(n1t,1));
UB = inf*ones(n1t,1);

[xu2,fo2,u2,yi2] = cplexmilp(c,-A1de,-yu,A1eq,b1eq,[],[],[],LB,UB,type,[],opt);
cys = [cys; [num2str(po), ') ', yi2.cplexstatusstring]]; 
if isempty(fo2)
    yi2;
    disp('Erro: Cplex não resolveu o problema inteiro.')
    continue
end

% fos e tempo
funp = [po fo1 fo2 tia yi1.time yi2.time];
funT = [funT; funp];

% avalia algumas porcentagens
if ~isempty(xu2)
    r = A2de*xu1;
    tok = (L.*W)'*xu2(1);
    pe = tok - ((r')*(l.*w));
    perda  = (pe*100)/tok;
    exce = (sum(r) - sum(d))*100;
    exc = exce/sum(d);
end

% informações das perdas
pep = [po perda exc  (sum(r) - sum(d))];
peT = [peT; pep];

% [W L]
% [w l d r]
disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = 'Demanda e produção';
tit = 'Itens. & $d_{i}$ & $pro_{i}$';
ProtoP = [(1:mwf)' d  r];
LatexTab(ProtoP,legdi,tit,[],1,[])
nam = num2str(po);
csvwrite(['CSVcom/',nam,'DemPro.csv'],ProtoP)

end


disp('========================================================' )
disp('RECORTAR AQUI')

disp(' ' )
disp(' ' )
disp('Problema, W, L, quantidade de itens, periodos') 
legdi = 'Informações do problema';
tit = 'Prob. & W & L & m & w & $\ell$ & d';
LatexTab(dork,legdi,tit,[],1,[0 0 0 0 1 0 1 0 1 0])
csvwrite('CSVcom/WLwld.csv',dork)


disp(' ' )
disp(' ' )
disp('Dimensão da matriz:') 
legdi = 'Dimensão da matriz do modelo';
tit = 'Prob. & nr1 & nv1 & nz1 & nr2 & nv2 & nz2';
disp('m da matriz, m max, n da matriz, n max') 
LatexTab(amax,legdi,tit,[ ],1,[])
csvwrite('CSVcom/mnModelo.csv',amax)


disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
tit = 'Prob. & f.o.i & T. mod. & T. Cplex ';
legdi = 'Dados da solução';
LatexTab(funT,legdi,tit,[ ],1,[])
csvwrite('CSVcom/FOs.csv',funT)


disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = 'Análise da produção';
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite('CSVcom/Perdas.csv',peT)

disp(cys)

diary off