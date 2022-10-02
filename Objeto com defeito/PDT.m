function [PT, feq, error] = PDT(w,ve,veO,rad,PO,PS)

[w,~]=sort(w);
mt = length(w);
error = 0;
vt = int16(full(ve));
vtO = int16(full(veO));
[lvetO,~,~]=find(vtO);
PT = [];

while ~isempty(lvetO)
    mlpO = length(lvetO);
    for ru =1:mlpO
        k = lvetO(ru);
        pd = zeros(mt,1);
        [lp,~,~]=find(PO(:,k));
        mlp = length(lp);
        for u = 1:mlp
            o = lp(u);
            col = PO(o,k);
            for j = 1:col
                certo = int16(vt(:,o));
                [lve,~,~] = find(certo>0);
                if ~isempty(lve)
                    if rad == 1
                    t = lve(1);
                    else
                    t = lve(randi(numel(lve)));
                    end
                else
                    error = 1;
                    [A,~,n] = unique(PT','rows');
                    [po,~] = size(A); 
                    feq = zeros(po,1);
                    PT = A';
                    for u = 1:po
                    [mmm,~,~]=find(n==u);
                    feq(u,1) = length(mmm);
                    end
                    return
                end
                pd = pd + PS(:,t,o);
                vt(t,o) = int16(vt(t,o)) - int16(1);
                vt(t,o) = int16(vt(t,o));
            end
        end
        vtO(k,1) = int16(vtO(k,1)) - int16(1);
        vtO(k,1) = int16(vtO(k,1));
        PT = [PT pd];
    end
[lvetO,~,~]=find(vtO);    
end

[A,~,n] = unique(PT','rows');
[po,~] = size(A); 
feq = zeros(po,1);

for u = 1:po
    [mmm,~,~]=find(n==u);
    feq(u,1) = length(mmm);
end

PT = A';

end



