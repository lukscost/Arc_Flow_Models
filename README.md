# Modelo de fluxo em arcos.

Códigos do MATLAB para resolver algumas versões do problema de corte de estoque bidimensional, guilhotinado, 2-estágios, com aparo (2D-PCE), utilizando o modelo de fluxo em arcos.

Para saber mais dos modelos ou para citação: (tese em andamento)  
Costa, L. L. S. (2022), Extensões do problema de corte de estoque bidimensional modelado como um problema de fluxo em arcos, Tese (doutorado), Universidade Estadual de Campinas, Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

# 3 - O modelo de fluxo em arcos

* 3.2 - Modelo de fluxo em arcos para o 2D-PCE.
Pasta: PCE com ou sem rotação
Arquivo principal: Bidimensional.m  

# 5 - Modelos propostos

* 5.1 - Modelo de fluxo em arcos para o 2D-PCE com múltiplos tipos de objetos  
Pasta: PCE com k-objetos  
Arquivo principal: kObjetos.m

* 5.2 - Modelo de fluxo em arcos para o 2D-PCE com múltiplos períodos
  * 5.2.2 - Modelo que permite estoque de itens  
  Pasta: PCE com periodos e estoque de itens  
  Arquivo principal: TPeriodos.m
  * 5.2.2.1 - Modelo que permite estoque de itens e com capacidade de estoque  
  Pasta: PCE com periodos e estoque de itens limite no estoque 
  Arquivo principal: TPeriodos.m

  * 5.2.3 - Modelo que permite estoque de faixas e itens  
Pasta: PCE com periodos e estoque de faixas 
Arquivo principal: TPeriodos.m

  * 5.2.4 - Modelo que permite atraso na produção de itens  
Pasta: PCE com periodos e atraso de itens 
Arquivo principal: TPeriodos.m

* 5.3 - Modelo de fluxo em arcos para o 2D-PCE com múltiplos tipos de objetos e múltiplos períodos  
  * 5.3.1 - Modelo que permite estoque de itens  
  Pasta: PCE com K-objetos e T-Períodos 
  Arquivo principal: KObjTPeri.m
  * 5.3.2 - Modelo que permite estoque de itens e com capacidade de estoque  
  Pasta: PCE com K-objetos e T-Períodos com limite de estoque
  Arquivo principal: KObjTPeri.m

* 5.4 - Modelo de fluxo em arcos para o 2D-PCE que permite rotação de itens  
  Pasta: PCE com ou sem rotação
  Arquivo principal: Bidimensional.m  

* 5.5 - Heurística de Decomposição do modelo de fluxo em arcos para o 2D-PCE
  * 5.5.x - Open Dimension 2D-PCE com rotação (Strip-Cutting Problem)
  Pasta: Decomposição  (com solução para o Strip-Cutting-Problem)
  * 5.5.1 - Com rotação, arquivo principal: BiComRotacao.m
  * 5.5.2 - Sem rotação, arquivo principal: BiSemRotacao.m
  Já apresentam o resultado da decomposição (Heurística) e do problema strip-cutting.  
 
* 5.6 - Adaptação do modelo de fluxo em arcos para construir padrões de corte bidimensional  
  * 5.6.1 - Modelo gerador de padrão de corte restrito  
  Pasta: Gerador de padrão com e sem rotação e restrito 
  Arquivo principal: Bidimensional.m
  Com e sem rotação, basta alterar o parâmetro rot: "hotrool = 1" (com rotação).  

  * 5.6.2 - Modelo gerador de padrão de corte irrestrito  
  Pasta: Gerador de padrão com e sem rotação e irrestrito 
  Arquivo principal: Bidimensional.m
  Com e sem rotação, basta alterar o parâmetro rot: "hotrool = 1" (com rotação).  

  * 5.6.3 - Adaptação do modelo de fluxo em arcos para construir padrões de corte bidimensional para tratar do caso com objeto com defeito.  
  Pasta: Objeto com defeito
  Arquivo principal: Layout












