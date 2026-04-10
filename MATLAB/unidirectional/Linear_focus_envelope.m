function linear_focus_envelope=Linear_focus_envelope(Akp,Alpha,x_vec,t,phishift)
%%%%%%%%%%%Hard-coded kp=0.0279 for Tp=12
kp=0.0279;
A=Akp/kp;
dx=x_vec(2)-x_vec(1);
Xmax=x_vec(end)-x_vec(1);
dk=2*pi/Xmax;
kmax=2*pi/dx/2;
N=length(x_vec);
k=linspace(0,kmax,floor(N/2)+1);
kw=0.004606;
g=9.81;
h=35;
t=1*t;
cp=sqrt(g/kp*tanh(kp*h));
cw=1/2*cp*(1+2*kp*h/sinh(2*kp*h));
for i=1:numel(k)
    if k(i)<kp
        S(i)=exp(-(k(i)-kp).^2/(2*kw.^2));
    else
        %        S(i)=exp(-(k(i)-kp).^2/(2*kw.^2))^(1/alpha);
        kw2=sqrt(kp^2/(2*log(10^Alpha)));
        S(i)=exp(-(k(i)-kp).^2/(2*kw2.^2));
    end
end
%         semilogy(k/kp,S);
%         ylim([1e-6,1e1]);
for i=1:numel(k)
    a(i)=A*S(i)*dk/(sum(S*dk));
end
w=sqrt(g.*k);
phi=-w*t+k.*t*cw;
% phi=-w*t;
F0=a.*numel(a).*exp(1i.*(phi+phishift));
F=[F0,conj(fliplr(F0(2:end-mod(N+1,2))))];
XX_linear=ifft(F);
XX_linear=fftshift(XX_linear);
%%%%%%%%%Compare with Gibbs
x_vec=linspace(0,dx*numel(XX_linear),numel(XX_linear));
x_vec=x_vec-x_vec(end)/2;
%                 Gibbs_envelope=A.*exp(-1/2.*kw^2.*x_vec.^2);
%                 hold on
%                 plot(Gibbs_envelope);
%                 plot(XX_linear);
linear_focus_envelope=(XX_linear);
end