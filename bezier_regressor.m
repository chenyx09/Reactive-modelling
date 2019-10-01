global bez_reg Bez_matr Bez_matr_der Bez_matr_dder N_bezier dB ddB
N_bezier=6;
dB=symbolic_dbezier(N_bezier);
ddB=symbolic_dbezier(N_bezier-1);
Nmin=8;
Nmax=200;
for i=Nmin:Nmax
    t=linspace(0,1,i);
    Bez_matr{i}=zeros(i,N_bezier+1);
    for j=1:i
        for k=0:N_bezier
            Bez_matr{i}(j,k+1)=nchoosek(N_bezier,k)*(1-t(j))^(N_bezier-k)*t(j)^k;
        end
        
    end
    bez_reg{i}=(Bez_matr{i}'*Bez_matr{i})\(Bez_matr{i}');
end


for i=Nmin:Nmax
    t=linspace(0,1,i);
    Bez_matr_der{i}=zeros(i,N_bezier);
    for j=1:i
        for k=0:N_bezier-1
            Bez_matr_der{i}(j,k+1)=nchoosek(N_bezier-1,k)*(1-t(j))^(N_bezier-1-k)*t(j)^k;
        end
    end
    
end

for i=Nmin:Nmax
    t=linspace(0,1,i);
    Bez_matr_dder{i}=zeros(i,N_bezier-1);
    for j=1:i
        for k=0:N_bezier-2
            Bez_matr_dder{i}(j,k+1)=nchoosek(N_bezier-2,k)*(1-t(j))^(N_bezier-2-k)*t(j)^k;
        end
    end
end