function P=gg_200921_numint_burst_1(Kini,b,g,beta,M,N,t,marg)
if strcmp(marg,'mature')
    M=1;
elseif strcmp(marg,'nascent')
    N=1;
end

l = 0:(M-1);
k = 0:(N-1);
u_ = exp(-2*pi*1i*l/M)-1;
v_ = exp(-2*pi*1i*k/N)-1;
[U,V] = ndgrid(u_,v_); u = U(:); v = V(:);

fun = @(x) INTFUNC(x,b,g,beta,u,v);
I___ = exp(Kini*integral(fun,0,t,'ArrayValued',true));

I___ = reshape(I___',[M,N])';
P=real(ifft2(I___))';
return

function F = INTFUNC(x,b,g,beta,U,V)
f = beta./(beta-g);
Ufun = f.*V.*exp(-g.*x) + (U-f.*V).*exp(-beta*x);
filt = ~isfinite(f);
if any(filt)
    Ufun(filt) = exp(-g*x).*(U+g.*V.*x);
end
Ufun = b.*Ufun;
F = Ufun./(1-Ufun);
return