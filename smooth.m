function out = smooth(in,n)
N=length(in);
d=(n-1)./2;
for i = 1:N
    if i >d && i<(N-d)
        out(i)=mean(in((i-d):(i+d)));
    elseif i<=d
        out(i)=mean(in(1:(i+d)));
    elseif i>=(N-d)
        out(i)=mean(in((i-d):N));
    end
end

