function f=heaviside(x)
    f=zeros(size(x));
    f(x==0)=1/2;
    f(x>0)=1;
end