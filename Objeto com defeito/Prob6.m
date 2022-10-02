% Costa, L. L. S. (2022), Extensões do problema de corte de estoque
% bidimensional modelado como um problema de fluxo em arcos,
% Tese (doutorado), Universidade Estadual de Campinas,
% Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

% Programa que trata do 2D-PCE, guilhotinado, 2-estágios e com aparo.
% Programa para problema com objeto com defeito

% Como a instância 6 obtivemos o seguinte resultado:
% {'6) Algum objeto residual não comporta nenhum item !' }
% Adaptamos o codigo para este caso, para considerar 3 objetos residuos

clc
clear all
close all
warning off 
format short

diary('Prob 6')
csvdir = 'CSV';
mkdir(csvdir)

opl
fou = [];
datis = []; 
amak = [];
dork = [];
cys = {};
Defyk = {};
peT = [];
ort = [];

disp('Ax >= d')

for jt= 6:6
gei = [];

disp(['Executando problema: ',num2str(jt),' ...']) 
% daddos
run('dados.m')

%  organização dos dados
[w,I]=sort(w); l=l(I); d = d(I); v = v(I);

% quantidade de itens, sem rotação e com 
mr = length(w);
% índices dos itens
ri = (1:mr)'; 
[mo,no] = size(WL);

disp(['Quantidade de itens: ',num2str(mr)])

PatK = [];
DemK = [];
rhs = [];
PatK1 = [];
DemK1 = [];
rhs1 = [];
vra = [];
ck = [];
yk = [];
amai = [0 0];

disp('Iniciando montagem da matriz do modelo ...')

try
tic;

for i=1:mo

if i == 3
    continue
end

ria = [];
lra = [];
wra = [];

for il =1:mr
  if (w(il) <= WL(i,1)) && (l(il) <= WL(i,2))
    ria = [ria; ri(il)];
    lra = [lra; l(il)];
    wra = [wra; w(il)];
  end
end

switch i
    case 1
        if isempty(ria)
         disp('Objeto gerado acima do defeito não comporta nenhum item !')
        Defyk  = [Defyk; [ num2str(jt),') Algum objeto residual não comporta nenhum item !']];
        end
         
    case 2
        if isempty(ria)
        disp('Objeto gerado abaixo do defeito não comporta nenhum item !')
        Defyk  = [Defyk; [ num2str(jt),') Algum objeto residual não comporta nenhum item !']];
        end

    case 3
        if isempty(ria)
        disp('Objeto gerado a esquerda do defeito não comporta nenhum item !')
        Defyk  = [Defyk; [ num2str(jt),') Algum objeto residual não comporta nenhum item !']];
        end

    case 4
        if isempty(ria)
        disp('Objeto gerado a direita do defeito não comporta nenhum item !')
        Defyk  = [Defyk; [ num2str(jt),') Algum objeto residual não comporta nenhum item !']];
        end
    otherwise
            disp('=============')
end

% criacao dos 'blocos duplos' correspondente aos dois objetos resultantens
% após a retirada do defeito
[DE, U, Tt, gt, itTO, ft, o4, idz, id1es, idzs , id2es] = MatMode(WL(i,2),lra,WL(i,1),wra,ria,mr);
[mt,nt] = size(PatK);
[mk,nk] = size(Tt);
[mt1,nt1] = size(PatK1);

% LAYOUT
switch i
    case 1
        PatK = Tt;
        DemK = DE;
        rhs = o4;
        c = zeros(nk,1);
        for ure = 1:mr
            Ire = find(itTO(id2es) == ure);
            c(id2es(Ire),1) = v(ure);
        end        
        ck = c;
    case 2
        zh1 = zeros(1,nt);
        zh1(1,1) = 1;
        zh2 = zeros(1,nk);
        zh2(1,1) = 1;
        vra = [vra zh1 zh2 0];
        yk = zeros(1,nk+nt+1);
        yk(1,nk+nt+1) = 1;
        PatK = [yk; zh1 zh2 -2; PatK zeros(mt,nk+1); zeros(mk,nt) Tt zeros(mk,1)];
        DemK = [DemK  DE zeros(mr,1)];
        rhs = [rhs; 0; o4];
        
        % vetor de custos do PI
        c = zeros(nk,1);
        for ure = 1:mr
            Ire = find(itTO(id2es) == ure);
            c(id2es(Ire),1) = v(ure);
        end
        ck = [ck; c; des(1,1)];
    case 3
        PatK1 = Tt;
        DemK1 = DE;
        rhs = [rhs; o4];
        % vetor de custos
        c = zeros(nk,1);
        for ure = 1:mr
            Ire = find(itTO(id2es) == ure);
            c(id2es(Ire),1) = v(ure);
        end
        ck = [ck; c];

    case 4
        zv1 = zeros(1,nt1);
        zv1(1,1) = 1;
        zv2 = zeros(1,nk);
        zv2(1,1) = 1;
        vra = [vra zv2 0];
        PatK1 = [zv2 -1; Tt zeros(mk,1)];
        DemK1 = [DemK1  DE zeros(mr,1)];
        rhs = [rhs; 0; o4];
        % vetor de custos do PI
        c = zeros(nk,1);
        for ure = 1:mr
            Ire = find(itTO(id2es) == ure);
            c(id2es(Ire),1) = v(ure);
        end
        ck = [ck; c; (des(1,2) - WL(3,1)*WL(3,2))];
        [mth,nth] = size(PatK);
	    [mtv,ntv] = size(PatK1);
        yv = zeros(1,nk+nt1+1);
        yv(1,nk+nt1+1) = 1;
        PatK = [PatK ([yv; zeros(mth-1,ntv)]); zeros(mtv,nth) PatK1];
        rhs = [1; rhs];
        DemK = [DemK  DemK1];
    otherwise
        PatK = [PatK zeros(mt,nk); zeros(mk,nt) Tt];
        DemK = [DemK  DE];
        rhs = [rhs; o4];
        % vetor de custos do PI
        c = zeros(nk,1);
        c(1) = WL(i,1)*WL(i,2);
        ck = [ck; c];
end
[mmax,nmax] = Ordem(WL(i,1),WL(i,2),l,w);
amai0 = [mmax nmax];
amai = amai + amai0;
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

% Para gerar os limitantes
[mt,nt] = size(PatK);
[mtd,ntd] = size(DemK);

datis0 = [jt mo  min(min(WL(:,1))) max(max(WL(:,1))) min(min(WL(:,2))) max(max(WL(:,2))) mr min(w) max(w) min(l) max(l) min(v) max(v)];
datis = [datis; datis0];

dens = (numel([PatK; DemK]) - nnz([PatK; DemK]))*100/numel([PatK; DemK]);
amak = [amak; jt (mt+mtd) amai(1,1) nt amai(1,2) dens];

% para checar a dimensão
[mck,nck] = size(ck);
[mt,nt] = size(PatK);
[mhs,nhs] = size(rhs);

[mdk,ndk] = size(DemK);
[mdd,ndd] = size(d);
[myk,nyk] = size([yk yv]);
 
% variáveis binárias
Iyk = find(PatK(1,:));
Izhv = find(vra);

LB = sparse(zeros(nt,1));
UB = inf*ones(nt,1);
UB(Iyk) = 1;
UB(Izhv) = 1;

tia = toc; 
disp(['Fez a matriz do modelo em ', num2str(tia)])
disp('Iniciando o CPLEX ...')

% tipo da solucao relaxada
opt = cplexoptimset('cplex');
opt.timelimit = 3*60; %8 min
% opt.display = 'on';

[x,fo,u,y] = cplexlp(-ck,[],[],PatK,rhs,LB,UB,[],opt);
y;
if isempty(fo)
    disp('Erro: Cplex não resolveu o problema relaxado.')
    continue
end

% quantidade de objetos usados
I = find(vra);
x(I);
y;

% tipo da solucao I = inteira, N = naturais
tipe = 'I'*ones(nt,1);
type = char(tipe');

[xu,foi,u,yi] = cplexmilp(-ck,[],[],PatK,rhs,[],[],[],LB,UB,type,[],opt);
% quantidade de objetos usados
cys = [cys; [ num2str(jt),') ', yi.cplexstatusstring]];
if isempty(foi)
    continue
end
obj = xu(I);
yi;
% disp('Solucao inteira')
foI = foi;
foi;

if xu(Iyk(1,1)) == 1 
    Defy = 'O DEFEITO FOI RETIRADO COM CORTES NA HORIZONTAL !';
    disp(Defy)
    disp(['O desperdício gerado para a retirada foi: ', num2str(des(1,1))])
    Defyk  = [Defyk; [ num2str(jt),') ', Defy]];
elseif xu(Iyk(1,2)) == 1
    Defy = 'O DEFEITO FOI RETIRADO COM CORTES NA VERTICAL !';
    disp(Defy)
    disp(['O desperdício gerado para a retirada foi: ', num2str(des(1,2))])
    Defyk  = [Defyk; [ num2str(jt),') ', Defy]];
else
    disp('O objeto com defeito não foi usado')
end

% avalia algumas porcentagens
if ~isempty(xu)
    r = DemK*xu;
    tok = W0*L0;
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

pep = [jt perda sum(r)];
peT = [peT; pep];
vapd = v'*r; 

disp(' ' )
disp(' ' )
disp('Quantidade de objetos usados:')

for ku = 1:3
    if obj(ku,1) ~=0
      dok = [ku WL(ku,1)  WL(ku,2) obj(ku,1)];
      gei = [gei; dok];
    else
      dok = [ku WL(ku,1)  WL(ku,2) 0];
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

disp(' ' )
disp('Itens, demanda-i, produção-i, i=1,2,3')
legdi = ['Demanda e produção; Colunas: 1) Itens; '...
    ' 2) Demanda do período;'...
    ' 3) Produção do período.'];
tit = 'Itens & $d_{i}$ & $Prod$';
ProtoP = [(1:mr)' w l v r];
LatexTab(ProtoP,legdi,tit,[],1,[ ])
nam = num2str(jt);
csvwrite([csvdir,'/',nam,'DemPro.csv'],ProtoP)

rtyo = ((r')*(l.*w));

fou0 = [jt vapd rtyo ((W0*L0)-rtyo+foi+vapd) (-foi-vapd) tia y.time];
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
tit = 'Prob. & k & min W & max W &  min L & max L  &  m & min w & max w & min l & max l &  min v & max v';
LatexTab(datis,legdi,tit,[],1,[ ])
csvwrite([csvdir,'/WLwld6.csv'],datis)


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
csvwrite([csvdir,'/mnModelo6.csv'],amak)

% fou
disp(' ' )
disp(' ' )
disp('total de objetos, objetos por periodo, fo, foi, tempo modelo, tempo cplex') 
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
csvwrite([csvdir,'/FOs6.csv'],fou)

disp(' ' )
disp(' ' )
disp('desperdicio da área total \% , exec produ \%, exec em itens') 
legdi = ['Análise da produção; Colunas: 1) Problema; '...
    ' 2) Porcentagem de desperdício da área total dos objetos;' ...
    ' 3) Porcentagem do excesso de produção da quantidade total de itens produzidos;' ...
    ' 4) Quantidade de itens produzidos a mais.'];
tit = 'Prob. & $(\%$ - desperdício$)$  & $(\%$ - exc. de prod.$)$  & exc. de itens';
LatexTab(peT,legdi,tit,[ ],1,[])
csvwrite([csvdir,'/Perdas6.csv'],peT)

cys
Defyk

diary off