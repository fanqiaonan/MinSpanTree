function [plv,const]=plv_hilbert(filtered1,filtered2)

%INPUT:filtered1,2 = 1 x n where n = number of samples
%signals should be filtered 
%OUTPUT: PLI VALUE
%LACHAUX "Measuring Phase Synchrony in Brain Signals" et al.1999

hilbert1=hilbert(filtered1);
hilbert2=hilbert(filtered2);

len=length(hilbert1);

phase1=angle(hilbert1);%������λֵ
phase2=angle(hilbert2);

unwrap1=unwrap(phase1);%����λ���ݽ��н����
unwrap2=unwrap(phase2);

plv=0;
plits=0; 

    for k=1:len
        const=unwrap1(k)-unwrap2(k);
        value=(exp(i*(unwrap1(k)-unwrap2(k))));
        plits=plits + value;
    end
 
    

    plv=abs(plits)/len;
