function []=Grafos(L, A, sca, capt, ref)

% full(A)

if isempty(sca)
sca = 0.25;
end
if isempty(capt)
    capt = 'Grafo do primeiro estagio';
end
if isempty(capt)
    ref = 'graf1';
end
    
%A = [2 0 4;
%     3 0 5;
%     3 5 10;
%     1 0 3;
%     1 3 6;
%     2 5 9
%     1 6 9];
% A = [ito, ga];

cor = {'blue', 'red', 'green', 'orange', 'magenta', 'cyan'};

[m,~] = size(A);

% define inicio
disp('\begin{figure}[H]')
disp('\centering')
disp(['\scalebox{',num2str(sca),'}{'])

% setup do grafo
disp('%\GraphInit[vstyle = Shade]')
disp('\tikzset{')
disp('LabelStyle/.style = {rectangle, rounded corners, draw, minimum width = 2em, font = \Large\bfseries},')
disp('VertexStyle/.append style = { inner sep=5pt,font = \Large\bfseries},')
disp('EdgeStyle/.append style = {->, bend left} }')

disp('\thispagestyle{empty}')

% grafo
disp('\begin{tikzpicture}')
disp('\SetGraphUnit{5}')

% inicio dos vertices
disp('\Vertex{0}')
disp('\EA(0){1}')
for i=1:(L-1)
    disp(['\EA(',num2str(i),'){',num2str(i+1),'}'])
end

% arcos
for j=1:m
% MUDA ENVERGADURA DO ARCO

if A(j,3) == 9999
    disp(['\tikzset{EdgeStyle/.append style = {bend left = 0', ...
        ', ultra thick, color=gray!50}}'])
    disp(['\Edge[label =', num2str(9999) ,'](', num2str(A(j,1)) ,')(',...
    num2str(A(j,2)) ,')'])
else
    disp(['\tikzset{EdgeStyle/.append style = {bend left =' , ...
        num2str(10 + 10*A(j,3)), ', ultra thick, color=', cor{A(j,3)},'!50}}'])
    disp(['\Edge[label =', num2str(A(j,3)) ,'](', num2str(A(j,1)) ,')(',...
        num2str(A(j,2)) ,')'])
end
    


end
 
% legenda  
disp('\end{tikzpicture}}')
disp(['\caption{',capt,'}'])
disp(['\label{',ref,'}'])
disp('\end{figure}')
















