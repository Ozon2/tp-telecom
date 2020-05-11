function [symboles] = A_PSKmod(bits, Nbits)
    symboles = zeros(1,Nbits/2);
    for k=1:Nbits/2
        bit1 = bits(2*k-1);
        bit2 = bits(2*k);
        if bit1 == 0 && bit2 == 0
            symboles(k)=-3;
        end
        if bit1 == 0 && bit2 == 1
            symboles(k)=-1;
        end
        if bit1 == 1 && bit2 == 1
            symboles(k)=1;                
        end
        if bit1 == 1 && bit2 == 0
            symboles(k)=3;
        end
    end
end

