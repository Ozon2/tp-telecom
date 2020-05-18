function [bits] = ASKdemod(symboles, Nbits)
    bits = zeros(1,Nbits);
    for k=1:Nbits/2
        if symboles(k) == -3
            bits(2*k-1)=0;
            bits(2*k)=0;
        end
        if symboles(k) == -1
            bits(2*k-1)=0;
            bits(2*k)=1;
        end
        if symboles(k) == 1
            bits(2*k-1)=1;
            bits(2*k)=1;               
        end
        if symboles(k) == 3
            bits(2*k-1)=1;
            bits(2*k)=0;
        end
    end
end
