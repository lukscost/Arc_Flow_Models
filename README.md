# Modelo de fluxo em arcos.

Códigos do MATLAB para resolver algumas versões do problema de corte de estoque bidimensional, guilhotinado, 2-estágios, com aparo, utilizando o modelo de fluxo em arcos.

Para saber mais dos modelos ou para citação: (tese em andamento)  
Costa, L. L. S. (202X), Extensões do problema de corte de estoque bidimensional modelado como um problema de fluxo em arcos, Tese (doutorado), Universidade Estadual de Campinas, Instituto de Matemática Estatística e Computação Científica, Campinas-SP.

# 3 - O modelo de fluxo em arcos

* 3.2 - Modelo de fluxo em arcos para o PCE bidimensional  
Pasta: Modelo de fluxo em arcos tradicional  
Arquivo principal: Bidimensional.m  

# 5 - Modelos propostos

* 5.1 - Modelo de fluxo em arcos para o PCE bidimensional com múltiplos tipos de objetos  
Pasta: Modelo com k tipos de objetos  
Arquivo principal: kObjetos.m

  * 5.1.1 - Adaptação do modelo de fluxo em arcos para o PCE bidimensional com múltiplos tipos de objetos para tratar do PCE bidimensional com alguns objetos defeituosos  
Em implementação...

* 5.2 - Modelo de fluxo em arcos para o PCE bidimensional com múltiplos períodos
  * 5.2.1 - Modelo que permite estoque de itens  
  Pasta: Modelo com mult. periodos e que permite estoque de itens  
  Arquivo principal: kPeriodos.m

  * 5.2.2 - Modelo que permite estoque de faixas e itens  
Pasta: Modelo com mult. periodos e que permite estoque de itens e faixas  
Arquivo principal: kPeriodos.m

  * 5.2.3 - Modelo que permite atraso na produção de itens  
Pasta: Modelo com mult. periodos e que permite atraso de itens  
Arquivo principal: kPeriodos.m

* 5.3 - Modelo de fluxo em arcos para o PCE bidimensional com múltiplos tipos de objetos e múltiplos períodos  
Em implementação ...

* 5.4 - Modelo de fluxo em arcos para o PCE bidimensional que permite rotação de itens  
Pasta: Modelo que permite rotação de itens  
Arquivo principal: Bidimensional.m  

* 5.5 - Decomposição do modelo de fluxo em arcos para o PCE bidimensional  
  * 5.5.1 - Open Dimension PCE bidimensional com rotação  
Pasta: Decomposição  
Permitindo rotação, arquivo principal: BiRotacao.m  
Sem permitir rotação, arquivo principal: Bidimensional.m  
Já apresentam o resultado da decomposição e do problema strip cutting.  
 
* 5.6 - Adaptação do modelo de fluxo em arcos para construir padrões de corte bidimensional  
  * 5.6.1 - Modelo gerador de padrão de corte restrito  
Pasta: Gerador de padrão  
Arquivo principal: Bidimensional.m  
Com e sem rotação, basta alterar o parâmetro rot: "rot = 1" (com rotação).  

  * 5.6.2 - Modelo gerador de padrão de corte irrestrito  
Não implementado pois, para os dados usados, alguns tipos de itens ocupam quase todo o objeto e por isso os padrões homogêneos só suportam um único item e não faz sentido restringir.











