function [xmx,xpx]=stftanalysis(x,w,N,H)
            %N:fft size,H:hop size
            b=size(w);
            M=b(1,2);
            hM1=floor((M+1)/2);
            hM2=floor(M/2);
            z1=zeros(1,hM1);
            z2=zeros(1,hM2);
            x=cat(2,z2,x);
            x=cat(2,x,z1);
            pin=hM1;
            b1=size(x);
            s1=b1(1,2);
            pend=s1-hM1;
            w=w/sum(w);
            xmx=[];
            xpx=[];
            while pin<=pend
                x1=x(1,pin-hM1+1:pin+hM2);
                [mx,px]=dftanalysis(x1,w,N);
                xmx=[xmx;mx];
                xpx=[xpx;px];
                pin=pin+H;
            end
end











