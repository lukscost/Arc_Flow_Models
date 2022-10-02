% criação dos dados para o problema com multiplos períodos

clear all
clc
diary('Data')

for j = 1:13
dar = ['Gcut\gcut',num2str(j),'d.m'];
run(dar)

a1 = 0.5;
a2 = 1.5;

r = a1 + (a2-a1).*rand(2,1);
d0 = [b (b.*r')] ;

dtab = [1 r'];
r';
d = ceil(d0);
dtab = [dtab' d'];

disp(['Dados do problema: ',num2str(j),':'])

disp(['Objetos: ',num2str(L),'x',num2str(W)])
disp('Itens: ')

legdi = ['Informação dos itens: 1) número do item, 2) largura e 3) comprimentos.'];
Co = num2str((1:3)');
gol = native2unicode(double(unicode2native(')'))*(ones(3,1)));
LatexTab([(1:length(w))' w l]',legdi,[],[Co gol],1)
disp('%==================')

disp('Demanda do problema: ')
legdi = ['Informação da demanda: 1) número do item, 2) demanda original, 3) e 4) demanda gerada multiplicando o valor da coluna (0).'];
Co = num2str((1:4)');
gol = native2unicode(double(unicode2native(')'))*(ones(4,1)));
LatexTab([[0 1:length(w)]; dtab],legdi,[],[Co gol],1)
disp('%==================')

nam = ['Peri',num2str(j),'.mat'];
save(nam,'L','l','W','w','d')

end

diary off



