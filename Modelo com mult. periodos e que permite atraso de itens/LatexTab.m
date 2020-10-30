function LatexTab(Ta,leg,tit,C,sc,co)

% leg = 'Tabela 1';
% tit = 'Dados & 1 & 2';
% X = randi(5,5);

[m,n]=size(Ta);

col = native2unicode(double(unicode2native('&'))*(ones(m,1)));
esp = native2unicode(double(unicode2native(' '))*(ones(m,1)));
ptvi = native2unicode(double(unicode2native(';'))*(ones(m,1)));
cop = native2unicode(double(unicode2native('['))*(ones(m,1)));
cclo = native2unicode(double(unicode2native(']'))*(ones(m,1)));
cf = native2unicode(double(unicode2native('$'))*(ones(m,1)));
nlin= native2unicode(double(unicode2native('\'))*(ones(m,1)));
XD = [];

if isempty(co)
    co = zeros(1,n);
end

i=1;
while i <= n
    if co(1,i) == 1
        XD = [XD cf cop num2str(Ta(:,i),'%4.2f') ptvi nlin esp  num2str(Ta(:,i+1),'%4.2f') cclo cf esp];    
        i = i+2;
    else
        XD = [XD num2str(Ta(:,i),'%4.2f') esp];    
        i = i+1;
    end
    if i <= n 
        XD = [XD  col esp];
    end
end


if ~isempty(C)
XD = [C esp col esp XD nlin nlin];
else
XD = [XD nlin nlin];    
end

disp('\begin{table}[!h]')
disp('\centering')
disp(['\scalebox{',num2str(sc),'}{ % tamanho da tabela'])

a =  [];
for j=1:n
    a =  [a 'r'];
end
a =  [a 'r'];

disp(['\begin{tabular}{', a ,'} % quantidade de colunas'])
% disp('\hline')

if ~isempty(tit)
    disp([tit,' \\'])
    disp('\toprule')
end

disp(XD)
disp('\bottomrule')
    

disp('\end{tabular}}')

if ~isempty(leg)
    disp(['\caption{', leg,' }'])
end

disp('\end{table}')

end