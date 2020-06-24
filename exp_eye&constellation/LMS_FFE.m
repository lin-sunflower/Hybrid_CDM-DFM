function [eq_out] = LMS_FFE(input,output,N,delta)
input=normalized_power(input);
output=normalized_power(output);

data1=input;
v=4;
for g=1:1:v
    data1=[data1 data1];
end
s=data1;

x=output(1:2.^v*length(input));
% figure;plot(x);hold on;plot(s)

out=x;re=s;
M=length(out);
x=out(1:N);
rxx=xcorr(x);
for i=1:N
    for j=1:N
        mrxx(i,j)=rxx(N-i+j);
    end
end

out1=zeros(1,M);
h=zeros(1,N);
for n=N:M
    x1=out(n:-1:n-N+1); 
    out1(n)=h*x1';
    e(n)=re(n-floor(N/2))-out1(n);
    h=h+delta.*e(n).*x1;
end
out_s=filter(h,1,output);
eq_out=out_s(floor(N/2)+1:end);
% out_s=filter(1,1,data);
% eq_out=out_s;
figure;plot(h)
end

