function [MST_LkMatrix] = Primf(CntMatrix)
    h = CntMatrix;
    l = size(CntMatrix,1);
    MST_LkMatrix = zeros(l,l);
    a = h';
    for i=1:l
        for j=1:l
             if a(i,j)==0
                   a(i,j)=inf;
             end
        end
    end
    k=1:l;
    listV(k)=0;
    listV(1)=1;
    e=1;
    while (e<l)
        min=inf;
         for i=1:l
            if listV(i)==1
                for j=1:l
                    if listV(j)==0
                       if min>a(i,j)
                            min=a(i,j);
                            b=a(i,j);
                            s=i;
                            d=j;
                        end
                    end
                end
            end
        end
        listV(d)=1;
        distance(e)=b;
        source(e)=s;
        destination(e)=d;
        e=e+1;
    end
    for g=1:e-1
        MST_LkMatrix(source(g),destination(g)) = distance(g);
        MST_LkMatrix(destination(g),source(g)) = distance(g);
    end
end