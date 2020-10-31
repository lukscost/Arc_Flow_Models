function lat(a,b,c)

[m,n]=size(a);

cf = native2unicode(double(unicode2native('$'))*(ones(m,1)));
col = native2unicode(double(unicode2native('&'))*(ones(m,1)));
esp = native2unicode(double(unicode2native(' '))*(ones(m,1)));
vig = native2unicode(double(unicode2native(','))*(ones(m,1)));
un = native2unicode(double(unicode2native('^'))*(ones(m,1)));
up = native2unicode(double(unicode2native('_'))*(ones(m,1)));
aco= native2unicode(double(unicode2native('{'))*(ones(m,1)));
fco= native2unicode(double(unicode2native('}'))*(ones(m,1)));
ax= native2unicode(double(unicode2native('x'))*(ones(m,1)));
ab= native2unicode(double(unicode2native('('))*(ones(m,1)));
fc = native2unicode(double(unicode2native(')'))*(ones(m,1)));

%Au = [cf ax up aco num2str(a,'%d') vig num2str(b,'%d') fco  un aco num2str(c,'%d') fco cf col];
Au = [cf ab num2str(a,'%d') vig num2str(b,'%d') fc cf col];
Bu = [];

for i=1:m
    Bu = [Bu Au(i,:)];
end


disp(Bu)
end