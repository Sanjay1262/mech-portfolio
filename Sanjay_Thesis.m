clear all
clc
%input paramter
P=2; %volume fraction index
Nu = 0.3;
rho0 = 1;
E0=10^9;
Em= 70*10^9;
Ec= 380*10^9;

mrho=2702;
crho=3800;
m = 1;
n = 1;
h = 1;
% t = 1;
mu=0.1; %Non-local parameter
c=2.5;
a=100*h;
b=a;
syms z
E=Em+(Ec-Em)*(0.5+z/h)^P;
D0=Ec*h^3/(12*(1-Nu^2));
K0=0;
K1=0;
k0=(K0*D0)/a^4;
k1=(K1*D0)/a^2;
syms x y 
lambda = (m*pi)/a;
Beta = (n*pi)/b;
L=1+mu*(lambda^2+Beta^2);
alphac=3.27*10^-6;
alpham=12.8*10^-6;
%%%%

%%%%%Displacement field
g=sin(c*z/h)*cos(c*z/h);  
Mu=(-c/h)*cos(c);
w=g+Mu*z;


h0=-h/2;
h1=-h/4;
h2=h/4;
h3=h/2;

%temperature gradients
delT1 = 0;
delT2 = 0;
delT3 = 0;

R1(z)= mrho +(crho-mrho)*((z-h0)/(h1-h0))^P;
R2= mrho +(crho-mrho)*1;
R3(z)= mrho +(crho-mrho)*((z-h3)/(h2-h3))^P;

Q111= Em/(1-Nu^2);
Q112= E/(1-Nu^2);
Q113= Ec/(1-Nu^2);
Q221= Em/(1-Nu^2);
Q222= E/(1-Nu^2);
Q223= Ec/(1-Nu^2);
Q661= Em/(2*(1+Nu));
Q662= E/(2*(1+Nu));
Q663= Ec/(2*(1+Nu));
Q121= Nu*Q111;
Q122= Nu*Q112;
Q123= Nu*Q113;
Q441= Em/(2*(1+Nu));
Q442= E/(2*(1+Nu));
Q443= Ec/(2*(1+Nu));

%%Matrix coeffecients
A11=(int(Q111,z,-h/2,-h/4) + int(Q112,z,-h/4,h/4) + int(Q113,z,h/4,h/2));
A22=(int(Q221,z,-h/2,-h/4) + int(Q222,z,-h/4,h/4) + int(Q223,z,h/4,h/2));
A66=(int(Q661,z,-h/2,-h/4) + int(Q662,z,-h/4,h/4) + int(Q663,z,h/4,h/2));
A12=(int(Q121,z,-h/2,-h/4) + int(Q122,z,-h/4,h/4) + int(Q123,z,h/4,h/2));

B11=(int(Q111*z,z,-h/2,-h/4) + int(Q112*z,z,-h/4,h/4) + int(Q113*z,z,h/4,h/2));
B22=(int(Q221*z,z,-h/2,-h/4) + int(Q222*z,z,-h/4,h/4) + int(Q223*z,z,h/4,h/2));
B66=(int(Q661*z,z,-h/2,-h/4) + int(Q662*z,z,-h/4,h/4) + int(Q663*z,z,h/4,h/2));
B12=(int(Q121*z,z,-h/2,-h/4) + int(Q122*z,z,-h/4,h/4) + int(Q123*z,z,h/4,h/2));

C11=(int(Q111*w,z,-h/2,-h/4) + int(Q112*w,z,-h/4,h/4) + int(Q113*w,z,h/4,h/2));
C22=(int(Q221*w,z,-h/2,-h/4) + int(Q222*w,z,-h/4,h/4) + int(Q223*w,z,h/4,h/2));
C66=(int(Q661*w,z,-h/2,-h/4) + int(Q662*w,z,-h/4,h/4) + int(Q663*w,z,h/4,h/2));
C12=(int(Q121*w,z,-h/2,-h/4) + int(Q122*w,z,-h/4,h/4) + int(Q123*w,z,h/4,h/2));

D11=(int(Q111*z^2,z,-h/2,-h/4) + int(Q112*z^2,z,-h/4,h/4) +int(Q113*z^2,z,h/4,h/2));
D22=(int(Q221*z^2,z,-h/2,-h/4) + int(Q222*z^2,z,-h/4,h/4) +int(Q223*z^2,z,h/4,h/2));
D66=(int(Q661*z^2,z,-h/2,-h/4) + int(Q662*z^2,z,-h/4,h/4) +int(Q663*z^2,z,h/4,h/2));
D12=(int(Q121*z^2,z,-h/2,-h/4) + int(Q122*z^2,z,-h/4,h/4) +int(Q123*z^2,z,h/4,h/2));

F11=(int(Q111*z*w,z,-h/2,-h/4) + int(Q112*z*w,z,-h/4,h/4) +int(Q113*z*w,z,h/4,h/2));
F22=(int(Q221*z*w,z,-h/2,-h/4) + int(Q222*z*w,z,-h/4,h/4) +int(Q223*z*w,z,h/4,h/2));
F66=(int(Q661*z*w,z,-h/2,-h/4) + int(Q662*z*w,z,-h/4,h/4) +int(Q663*z*w,z,h/4,h/2));
F12=(int(Q121*z*w,z,-h/2,-h/4) + int(Q122*z*w,z,-h/4,h/4) +int(Q123*z*w,z,h/4,h/2));

H11=(int(Q111*w^2,z,-h/2,-h/4) + int(Q112*w^2,z,-h/4,h/4) +int(Q113*w^2,z,h/4,h/2));
H22=(int(Q221*w^2,z,-h/2,-h/4) + int(Q222*w^2,z,-h/4,h/4) +int(Q223*w^2,z,h/4,h/2));
H66=(int(Q661*w^2,z,-h/2,-h/4) + int(Q662*w^2,z,-h/4,h/4) +int(Q663*w^2,z,h/4,h/2));
H12=(int(Q121*w^2,z,-h/2,-h/4) + int(Q122*w^2,z,-h/4,h/4) +int(Q123*w^2,z,h/4,h/2));

J44=(int(Q441*diff(w)^2,z,-h/2,-h/4) + int(Q442*diff(w)^2,z,-h/4,h/4) +int(Q443*diff(w)^2,z,h/4,h/2));
J66=(int(Q661*diff(w)^2,z,-h/2,-h/4) + int(Q662*diff(w)^2,z,-h/4,h/4) +int(Q663*diff(w)^2,z,h/4,h/2));

At11=(int((alphac*Em*delT1)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2)/(1-Nu),z,-h/4,h/4) +int((alphac*Em*delT3)/(1-Nu),z,h/4,h/2));
At22=(int((alphac*Em*delT1)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2)/(1-Nu),z,-h/4,h/4) +int((alphac*Em*delT3)/(1-Nu),z,h/4,h/2));

Bt11=(int((alphac*Em*delT1*z)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z)/(1-Nu),z,h/4,h/2));
Bt22=(int((alphac*Em*delT1*z)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z)/(1-Nu),z,h/4,h/2));

Dt11=(int((alphac*Em*delT1*z^2)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z^2)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z^2)/(1-Nu),z,h/4,h/2));
Dt22=(int((alphac*Em*delT1*z^2)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z^2)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z^2)/(1-Nu),z,h/4,h/2));

Ct11=(int((alphac*Em*delT1*w)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*w)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*w)/(1-Nu),z,h/4,h/2));
Ct22=(int((alphac*Em*delT1*w)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*w)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*w)/(1-Nu),z,h/4,h/2));

Ft11=(int((alphac*Em*delT1*z*w)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z*w)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z*w)/(1-Nu),z,h/4,h/2));
Ft22=(int((alphac*Em*delT1*z*w)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*z*w)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*z*w)/(1-Nu),z,h/4,h/2));

Ht11=(int((alphac*Em*delT1*w^2)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*w^2)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*w^2)/(1-Nu),z,h/4,h/2));
Ht22=(int((alphac*Em*delT1*w^2)/(1-Nu),z,-h/2,-h/4)+int((alphac*Ec*delT2*w^2)/(1-Nu),z,-h/4,h/4)+int((alphac*Em*delT3*w^2)/(1-Nu),z,h/4,h/2));

I0=(int(R1,z,-h/2,-h/4) + int(R2,z,-h/4,h/4) + int(R3,z,h/4,h/2));
I1=(int(z*R1,z,-h/2,-h/4) + int(z*R2,z,-h/4,h/4) + int(z*R3,z,h/4,h/2));
I2=(int(z^2*R1,z,-h/2,-h/4) + int(z^2*R2,z,-h/4,h/4) + int(z^2*R3,z,h/4,h/2));
I3=(int(w*R1,z,-h/2,-h/4) + int(w*R2,z,-h/4,h/4) + int(w*R3,z,h/4,h/2));
I4=(int(z*w*R1,z,-h/2,-h/4) + int(z*w*R2,z,-h/4,h/4) + int(z*w*R3,z,h/4,h/2));
I5=(int(w^2*R1,z,-h/2,-h/4) + int(w^2*R2,z,-h/4,h/4) + int(w^2*R3,z,h/4,h/2));
%%%%

[ty1,ty2,ty3,ty4]=deal(0);
for m=1:1
    for n=1:1
        if bitget(m,1)&&bitget(n,1)
            lambda=m*pi/a;
            Beta=n*pi/b;
            K=zeros(5,5);
            K(1,1) =K(1,1) +(A11+ At11)*lambda^2 + (A66+At22)*Beta^2;
            K(1,2) =lambda*Beta*(A12+A66);
            K(1,3) =-(B11+Bt11)*lambda^3 - lambda*(Beta)^2*(B12+2*B66+Bt22);
            K(1,4) =(C11+Ct11)*lambda^2 +(C66+Ct22)*Beta^2;
            K(1,5) =lambda*Beta*(C12+C66);
            K(2,2) =(A66+ At11)*lambda^2 + (A22+At22)*Beta^2;
            K(2,5) =(C66+Ct11)*lambda^2 +(C22+Ct22)*Beta^2;
            K(2,3) =-(B22+Bt22)*Beta^3 - Beta*(lambda)^2*(B12+2*B66+Bt11);
            K(2,4) =K(1,5);
            K(3,3) =(D11+Dt11)*lambda^4+(2*D12+4*D66+Dt11+Dt22)*(lambda^2)*(Beta^2)+(D22+Dt22)*Beta^4+L*(At11*(lambda)^2+At22*(Beta)^2);
            K(3,4) =-(F11+Ft11)*lambda^3-lambda*(Beta)^2*(F12+2*F66+Ft22);
            K(3,5) =-(F22+Ft22)*Beta^3-Beta*(lambda)^2*(F12+2*F66+Ft11);
            K(4,4) =J44+(H11+Ht11)*lambda^2+(H66+Ht22)*Beta^2;
            K(4,5) =lambda*Beta*(H12+H66);
            K(5,5) =J44+(H66+Ht11)*lambda^2+(H22+Ht22)*Beta^2;
            q0=1;
            Qmn=q0;
            load=[0 0 Qmn 0 0]';
            u=K\load;
            ry1=vpa(u(1)*cos(m*pi*x/a)*sin(n*pi*y/a));
            ry2=vpa(u(2)*sin(m*pi*x/a)*cos(n*pi*y/b));
            ry3=vpa(u(3)*sin(m*pi*x/a)*sin(n*pi*y/b));
            ry4=vpa(u(4)*sin(m*pi*x/a)*sin(n*pi*y/b));
            ty1=ty1+ry1;
            ty2=ty2+ry2;
            ty3=ty3+ry3;
            ty4=ty4+ry4;
        end
    end
end
U_DEF=ty1;
V_DEF=ty2;
Wb_DEF=ty3;
Ws_DEF=ty4;
w_def=ty3+ty4;
Non_WDef=double(10*Ec*h^3*subs(w_def,{x,y},{a/2,b/2})/(a^4*q0));
m=1; n=1;
K= zeros(5,5);
M= zeros(5,5);
K(1,1) = (A11+ At11)*lambda^2 + (A66+At22)*Beta^2;
K(1,2) = lambda*Beta*(A12+A66);
K(1,3) = -(B11+Bt11)*lambda^3 - lambda*(Beta)^2*(B12+2*B66+Bt22);
K(1,4) = (C11+ Ct11)*lambda^2 + (C66 + Ct22)*Beta^2;
K(1,5) = lambda*Beta*(C12+C66);
K(2,2) = (A66+ At11)*lambda^2 + (A22+At22)*Beta^2;
K(2,3) = -(B22+Bt22)*Beta^3 - Beta*(lambda)^2*(B12+2*B66+Bt11);
K(2,4) = K(1,5);
K(2,5) = (C66+Ct11)*lambda^2 +(C22+Ct22)*Beta^2;
K(3,3) = (D11+Dt11)*lambda^4+(2*D12+4*D66+Dt11+Dt22)*(lambda^2)*(Beta^2)+(D22+Dt22)*Beta^4+k0+k1*L*(At11*(lambda)^2+At22*(Beta)^2);
K(3,4) = -(F11+Ft11)*lambda^3-lambda*(Beta)^2*(F12+2*F66+Ft22);
K(3,5) = -(F22+Ft22)*Beta^3-Beta*(lambda)^2*(F12+2*F66+Ft11);
K(4,4) = J44+(H11+Ht11)*lambda^2+(H66+Ht22)*Beta^2;
K(4,5) = lambda*Beta*(H12+H66);
K(5,5) = J44+(H66+Ht11)*lambda^2+(H22+Ht22)*Beta^2;
M(1,1) = L*I0;
M(1,3) = -lambda*L*I1;  
M(1,4) = L*I3;
M(1,5) = L*I3;
M(2,2) = L*I0;
M(2,3) = -Beta*L*I1;
M(3,3) = (L*I0)-(L*I2*(lambda^2+Beta^2));
M(3,4) = -lambda*L*I4;
M(3,5) = -Beta*L*I4;
M(4,4) = L*I5;
M(5,5) = L*I5;
M(1,2)=0; M(1,4)=0; M(1,5)=0; M(4,5)=0;
omn = eig(K,M);
frqm = abs(sqrt(omn));
frq1 = (frqm*(a^2/h)*sqrt(rho0/E0));
frq=(10*frqm*sqrt(mrho/Em));
disp(min(frq1));

 c = 2.5;
% % g = (c*z)/(h*(((c^2*z^2))/h^2)+1);
% % Mu = -(8*c*h^2)/(c^3*h^3+8*h^3);
 g=sin(c*z/h)*cos(c*z/h);
 Mu=(-c/h)*cos(c);
 f = g + Mu*z;
 
 phix = diff(diff(Ws_DEF, x));
 phiy = diff(diff(Wb_DEF, y));
 
 u = U_DEF - z*diff(Wb_DEF, x) + f*phix;
 v = V_DEF - z*diff(Wb_DEF, y) + f*phiy;
 w = Wb_DEF;
 
 epsx=diff(U_DEF,x)-z*diff(diff(Wb_DEF,x),x)+[g + Mu*z]*(phix);
 epsy=diff(V_DEF,y)-z*diff(diff(Wb_DEF,y),y)+[g + Mu*z]*(phiy);
 gamaxy=diff(U_DEF,y)+diff(V_DEF,x)-2*z*diff(diff(Wb_DEF,x),y)+[g + Mu*z]*(diff(phix,y)+diff(phiy,x));
 gamaxz= (diff(g, z) + Mu)*diff(Ws_DEF, x);
 Gammayz = (diff(g, z) + Mu)*diff(Wb_DEF, y);
 %Z = z_values/h;
% % 
% % % Plot
% % plot(Gammayz_values_double, Z);
% % xlabel('Gamma yz');
% % ylabel('Z/h');
% % epsx=diff(U_DEF,x)-Z*diff(diff(Wb_DEF,x))-4*Z^3/(3*h^2)*diff(diff(Ws_DEF,x));
% % epsy=diff(V_DEF,y)-Z*diff(diff(Wb_DEF,y))-4*Z^3/(3*h^2)*diff(diff(Ws_DEF,y));
% % gamaxy=diff(u(2),x)+diff(u(1),y)-2*Z*diff((diff(Wb_DEF,x)),y)-8*Z^3/(3*h^2)*diff((diff(Ws_DEF,x)),y);
% % gamaxz=(1-4*Z^2/h^2)*diff(Ws_DEF,x);
% % gamayz=(1-4*Z^2/h^2)*diff(Wb_DEF,y);
 strain=[epsx epsy gamaxy gamaxz Gammayz]';
 E(z)=Em+(Ec-Em)*(0.5+z/h)^P;
 coef=E(z)/(1-Nu)*[1 Nu 0 0 0; Nu 1 0 0 0; 0 0 (1-Nu)/2 0 0; 0 0 0 (1-Nu)/2 0; 0 0 0 0 (1-Nu)/2];
 sigma=coef*strain;
 sigmazz=subs(sigma,{x,y},{a/2,b/2});
 sigmayz=h/(q0*a)*sigmazz(4 );
 tt=-0.5:0.01:0.5;
 tt = h/6; % Set the z value for calculation
 kk1 = double(subs(sigmayz, z, tt));
 
 
 z_values = linspace(-h/2, h/2); 
 
 kk = subs(sigmayz, z, z_values);
 
 % Plot the results
 plot(kk, z_values); 
 xlabel('Sigma x'); 
 ylabel('Z/h');
 title('Stress Sigma x vs. Z/h'); 
 grid on;