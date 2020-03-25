function M=MCC_coeff(C)
[~,n]=size(C);
nomin=0;
temp=0;
denomin_1=0;
denomin_2=0;

for k=1:n
    for l=1:n
        for m=1:n
            nomin=nomin+C(k,k)*C(m,l)-C(l,k)*C(k,m);
        end
    end
end

for k=1:n
    for l=1:n
        for f=1:n
            for g=1:n
                if f~=k
                   temp=temp+C(g,f);                       
                end
            end 
        end
        denomin_1=denomin_1+C(l,k)*temp;
        temp=0;
    end
end



for k=1:n
    for l=1:n
        for f=1:n
            for g=1:n
                if f~=k
                    temp=temp+C(f,g);
                end
            end
        end
        denomin_2=denomin_2+C(k,l)*temp;
        temp=0;
    end
end

M=nomin/sqrt(denomin_1*denomin_2);

if nomin==0 
    M=0;
end

end
            
        