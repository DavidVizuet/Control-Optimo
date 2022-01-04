function y = SEIR(a,b,c,d,e,g,S0,E0,I0,R0,T)

test = -1;

delta = 0.001;
M = 1000;
t=linspace(0,T,M+1);
h=T/M;
h2 = h/2;

S=zeros(1,M+1);
E=zeros(1,M+1);
I=zeros(1,M+1);
R=zeros(1,M+1);
N=zeros(1,M+1);
S(1)=S0;
E(1)=E0;
I(1)=I0;
R(1)=R0;
N(1)=S0+E0+I0+R0;

for i = 1:M
        
        n11 = b*N(i) - d*S(i) - c*S(i)*I(i);
        n12 = c*S(i)*I(i) - (e+d)*E(i);
        n13 = e*E(i) - (g+a+d)*I(i);
        n14 = (b-d)*N(i) - a*I(i);
        
        n21 = b*(N(i)+h2*n14) - d*(S(i)+h2*n11) - c*(S(i)+h2*n11)*(I(i)+h2*n13);
        n22 = c*(S(i)+h2*n11)*(I(i)+h2*n13) - (e+d)*(E(i)+h2*n12);
        n23 = e*(E(i)+h2*n12) - (g+a+d)*(I(i)+h2*n13);
        n24 = (b-d)*(N(i)+h2*n14) - a*(I(i)+h2*n13);
        
        n31 = b*(N(i)+h2*n24) - d*(S(i)+h2*n21) - c*(S(i)+h2*n21)*(I(i)+h2*n23);
        n32 = c*(S(i)+h2*n21)*(I(i)+h2*n23) - (e+d)*(E(i)+h2*n22);
        n33 = e*(E(i)+h2*n22) - (g+a+d)*(I(i)+h2*n23);
        n34 = (b-d)*(N(i)+h2*n24) - a*(I(i)+h2*n23);
        
        n41 = b*(N(i)+h*n34) - d*(S(i)+h*n31) - c*(S(i)+h*n31)*(I(i)+h*n33);
        n42 = c*(S(i)+h*n31)*(I(i)+h*n33) - (e+d)*(E(i)+h*n32);
        n43 = e*(E(i)+h*n32) - (g+a+d)*(I(i)+h*n33);
        n44 = (b-d)*(N(i)+h*n34) - a*(I(i)+h*n33);
        
        S(i+1) = S(i) + (h/6)*(n11 + 2*n21 + 2*n31 + n41);
        E(i+1) = E(i) + (h/6)*(n12 + 2*n22 + 2*n32 + n42);
        I(i+1) = I(i) + (h/6)*(n13 + 2*n23 + 2*n33 + n43);
        N(i+1) = N(i) + (h/6)*(n14 + 2*n24 + 2*n34 + n44);
end

for i=1:M
    n1 = g*I(i) - d*R(i); 
    n2 = g*0.5*(I(i)+I(i+1)) - d*(R(i)+h2*n1);
    n3 = g*0.5*(I(i)+I(i+1)) - d*(R(i)+h2*n2); 
    n4 = g*I(i+1) - d*(R(i)+h*n3);
    R(i+1) = R(i) + (h/6)*(n1 + 2*n2 + 2*n3 + n4);
end

y(1,:) = t;
y(2,:) = S;
y(3,:) = E;
y(4,:) = I;
y(5,:) = R;
y(6,:) = N;

plot(y(1,:),y(6,:),'LineWidth',3) 
xlabel('tiempo (en años)')
title("Población total")

%plot(y(1,:),y(3,:)) 
%title("Población Total")

%plot(y(1,:),y(4,:)) 
%title("Infectados")

%plot(y(1,:),y(5,:)) 
%title("Recuperados")

%plot(y(1,:),y(6,:)) 
%title("Todos")


%subplot(3,2,1),plot(y(1,:),y(2,:),'LineWidth',2),title("Susceptibles");
%subplot(3,2,2),plot(y(1,:),y(3,:),'LineWidth',2),title("Expuestos");
%subplot(3,2,3),plot(y(1,:),y(4,:),'LineWidth',2),title("Infectados");
%subplot(3,2,4),plot(y(1,:),y(5,:),'LineWidth',2),title("Recuperados");
%subplot(3,2,5),plot(y(1,:),y(6,:),'LineWidth',2),title("Población TOTAL");
%subplot(3,2,6),plot(y(1,:),y(7,:),'LineWidth',2),title("Control");

