function y = Control(b,d,c,e,g,a,S0,E0,I0,R0,A,T)

test = -1;

delta = 0.001;
M = 1000;
t=linspace(0,T,M+1);
h=T/M;
h2 = h/2;

%Variables en el modelo SEIR 
Sn=zeros(1,M+1);
En=zeros(1,M+1);
In=zeros(1,M+1);
Rn=zeros(1,M+1);
Nn=zeros(1,M+1);
Sn(1)=S0;
En(1)=E0;
In(1)=I0;
Rn(1)=R0;
Nn(1)=S0+E0+I0+R0;

%Variables en el modelo SEIR con control
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


%Variables adjuntas del sistema adjunto
lambda1=zeros(1,M+1);
lambda2=zeros(1,M+1);
lambda3=zeros(1,M+1);
lambda4=zeros(1,M+1);

u=zeros(1,M+1);

while(test < 0)
    
    oldu = u;
    oldS = S;
    oldE = E;
    oldI = I;
    oldN = N;
    oldlambda1 = lambda1;
    oldlambda2 = lambda2;
    oldlambda3 = lambda3;
    oldlambda4 = lambda4;
    
    for i = 1:M
        % Implementacion del metodo de Runge-Kutta para el modelo SEIR 
        n11 = b*Nn(i) - d*Sn(i) - c*Sn(i)*In(i);
        n12 = c*Sn(i)*In(i) - (e+d)*En(i);
        n13 = e*En(i) - (g+a+d)*In(i);
        n14 = (b-d)*Nn(i) - a*In(i);
        
        n21 = b*(Nn(i)+h2*n14) - d*(Sn(i)+h2*n11) - c*(Sn(i)+h2*n11)*(In(i)+h2*n13);
        n22 = c*(Sn(i)+h2*n11)*(In(i)+h2*n13) - (e+d)*(En(i)+h2*n12);
        n23 = e*(En(i)+h2*n12) - (g+a+d)*(In(i)+h2*n13);
        n24 = (b-d)*(Nn(i)+h2*n14) - a*(In(i)+h2*n13);
        
        n31 = b*(Nn(i)+h2*n24) - d*(Sn(i)+h2*n21) - c*(Sn(i)+h2*n21)*(In(i)+h2*n23);
        n32 = c*(Sn(i)+h2*n21)*(In(i)+h2*n23) - (e+d)*(En(i)+h2*n22);
        n33 = e*(En(i)+h2*n22) - (g+a+d)*(In(i)+h2*n23);
        n34 = (b-d)*(Nn(i)+h2*n24) - a*(In(i)+h2*n23);
        
        n41 = b*(Nn(i)+h*n34) - d*(Sn(i)+h*n31) - c*(Sn(i)+h*n31)*(In(i)+h*n33);
        n42 = c*(Sn(i)+h*n31)*(In(i)+h*n33) - (e+d)*(En(i)+h*n32);
        n43 = e*(En(i)+h*n32) - (g+a+d)*(In(i)+h*n33);
        n44 = (b-d)*(Nn(i)+h*n34) - a*(In(i)+h*n33);
        
        Sn(i+1) = Sn(i) + (h/6)*(n11 + 2*n21 + 2*n31 + n41);
        En(i+1) = En(i) + (h/6)*(n12 + 2*n22 + 2*n32 + n42);
        In(i+1) = In(i) + (h/6)*(n13 + 2*n23 + 2*n33 + n43);
        Nn(i+1) = Nn(i) + (h/6)*(n14 + 2*n24 + 2*n34 + n44);
        
        
        
        
        %Implementacion del metodo de Runge Kutta al resolver las variables de estado hacia adelante en el tiempo
        m11 = b*N(i) - d*S(i) - c*S(i)*I(i) - u(i)*S(i);
        m12 = c*S(i)*I(i) - (e+d)*E(i);
        m13 = e*E(i) - (g+a+d)*I(i);
        m14 = (b-d)*N(i) - a*I(i);
        
        m21 = b*(N(i)+h2*m14) - d*(S(i)+h2*m11) - c*(S(i)+h2*m11)*(I(i)+h2*m13) - 0.5*(u(i)+u(i+1))*(S(i)+h2*m11);
        m22 = c*(S(i)+h2*m11)*(I(i)+h2*m13) - (e+d)*(E(i)+h2*m12);
        m23 = e*(E(i)+h2*m12) - (g+a+d)*(I(i)+h2*m13);
        m24 = (b-d)*(N(i)+h2*m14) - a*(I(i)+h2*m13);
        
        m31 = b*(N(i)+h2*m24) - d*(S(i)+h2*m21) - c*(S(i)+h2*m21)*(I(i)+h2*m23) - 0.5*(u(i)+u(i+1))*(S(i)+h2*m21);
        m32 = c*(S(i)+h2*m21)*(I(i)+h2*m23) - (e+d)*(E(i)+h2*m22);
        m33 = e*(E(i)+h2*m22) - (g+a+d)*(I(i)+h2*m23);
        m34 = (b-d)*(N(i)+h2*m24) - a*(I(i)+h2*m23);
        
        m41 = b*(N(i)+h*m34) - d*(S(i)+h*m31) - c*(S(i)+h*m31)*(I(i)+h*m33) - u(i+1)*(S(i)+h*m31);
        m42 = c*(S(i)+h*m31)*(I(i)+h*m33) - (e+d)*(E(i)+h*m32);
        m43 = e*(E(i)+h*m32) - (g+a+d)*(I(i)+h*m33);
        m44 = (b-d)*(N(i)+h*m34) - a*(I(i)+h*m33);
        
        S(i+1) = S(i) + (h/6)*(m11 + 2*m21 + 2*m31 + m41);
        E(i+1) = E(i) + (h/6)*(m12 + 2*m22 + 2*m32 + m42);
        I(i+1) = I(i) + (h/6)*(m13 + 2*m23 + 2*m33 + m43);
        N(i+1) = N(i) + (h/6)*(m14 + 2*m24 + 2*m34 + m44);
    end

    
    %Implementacion del metodo de Runge Kutta resolviendo el sistema adjunto (sistema para las lambdas) hacia atras en el tiempo
    for i = 1:M
        j = M + 2 - i;
        m11 = lambda1(j)*(d + c*I(j) + u(j)) - c*lambda2(j)*I(j);
        m12 = lambda2(j)*(e + d) - lambda3(j)*e;
        m13 = -A + (lambda1(j) - lambda2(j))*c*S(j) + lambda3(j)*(g+a+d) + lambda4(j)*a;
        m14 = -lambda1(j)*b - lambda4(j)*(b-d);
        
        m21 = (lambda1(j)-h2*m11)*(d + c*0.5*(I(j) + I(j-1)) + 0.5*(u(j) + u(j-1))) - c*(lambda2(j)-h2*m12)*0.5*(I(j) + I(j-1));
        m22 = (lambda2(j)-h2*m12)*(e + d) - (lambda3(j)-h2*m13)*e;
        m23 = -A + ((lambda1(j)-h2*m11) - (lambda2(j)-h2*m12))*c*0.5*(S(j) + S(j-1)) + (lambda3(j)-h2*m13)*(g+a+d) + (lambda4(j)-h2*m14)*a;
        m24 = -(lambda1(j)-h2*m11)*b - (lambda4(j)-h2*m14)*(b-d);
        
        m31 = (lambda1(j)-h2*m21)*(d + c*0.5*(I(j) + I(j-1)) + 0.5*(u(j) + u(j-1))) - c*(lambda2(j)-h2*m22)*0.5*(I(j) + I(j-1));
        m32 = (lambda2(j)-h2*m22)*(e + d) - (lambda3(j)-h2*m23)*e;
        m33 = -A + ((lambda1(j)-h2*m21) - (lambda2(j)-h2*m22))*c*0.5*(S(j) + S(j-1)) + (lambda3(j)-h2*m23)*(g+a+d) + (lambda4(j)-h2*m24)*a;
        m34 = -(lambda1(j)-h2*m21)*b - (lambda4(j)-h2*m24)*(b-d);
        
        m41 = (lambda1(j)-h*m31)*(d + c*I(j-1) + u(j-1)) - c*(lambda2(j)-h*m32)*I(j-1);
        m42 = (lambda2(j)-h*m32)*(e + d) - (lambda3(j)-h*m33)*e;
        m43 = -A + ((lambda1(j)-h*m31) - (lambda2(j)-h*m32))*c*S(j-1) + (lambda3(j)-h*m33)*(g+a+d) + (lambda4(j)-h*m34)*a;
        m44 = -(lambda1(j)-h*m31)*b - (lambda4(j)-h*m34)*(b-d);
        
        lambda1(j-1) = lambda1(j) - (h/6)*(m11 + 2*m21 + 2*m31 + m41);
        lambda2(j-1) = lambda2(j) - (h/6)*(m12 + 2*m22 + 2*m32 + m42);
        lambda3(j-1) = lambda3(j) - (h/6)*(m13 + 2*m23 + 2*m33 + m43);
        lambda4(j-1) = lambda4(j) - (h/6)*(m14 + 2*m24 + 2*m34 + m44);
    end
        %Representacion de las variables de control u(t) usando los nuevos valores para lambda
        temp=(S.*lambda1)./2;
        u1 = min(0.9,max(0,temp));
        u = 0.5*(u1 + oldu);
       
       
    %Parametros para estimar la convergencia de cada variable. 
    temp1 = delta*sum(abs(u)) - sum(abs(oldu - u));
    temp2 = delta*sum(abs(S)) - sum(abs(oldS - S));
    temp3 = delta*sum(abs(E)) - sum(abs(oldE - E));
    temp4 = delta*sum(abs(I)) - sum(abs(oldI - I));
    temp5 = delta*sum(abs(N)) - sum(abs(oldN - N));
    temp6 = delta*sum(abs(lambda1)) - sum(abs(oldlambda1 - lambda1));
    temp7 = delta*sum(abs(lambda2)) - sum(abs(oldlambda2 - lambda2));
    temp8 = delta*sum(abs(lambda3)) - sum(abs(oldlambda3 - lambda3));
    temp9 = delta*sum(abs(lambda4)) - sum(abs(oldlambda4 - lambda4));
    test = min(temp1, min(temp2, min(temp3, min(temp4, min(temp5, min(temp6, min(temp7, min(temp8, temp9))))))));
end

%Con los valores obtenidos de la variable u(t) se implementa el metodo de Runge Kutta para resolver la ecuacion de los recuperados
for i=1:M
    n1 = g*In(i) - d*Rn(i); 
    n2 = g*0.5*(In(i)+In(i+1)) - d*(Rn(i)+h2*n1);
    n3 = g*0.5*(In(i)+In(i+1)) - d*(Rn(i)+h2*n2); 
    n4 = g*In(i+1) - d*(Rn(i)+h*n3);
    Rn(i+1) = Rn(i) + (h/6)*(n1 + 2*n2 + 2*n3 + n4);
    
    m1 = g*I(i) - d*R(i) + u(i)*S(i);
    m2 = g*0.5*(I(i)+I(i+1)) - d*(R(i)+h2*m1) + 0.5*(u(i)+u(i+1))*0.5*(S(i)+S(i+1));
    m3 = g*0.5*(I(i)+I(i+1)) - d*(R(i)+h2*m2) + 0.5*(u(i)+u(i+1))*0.5*(S(i)+S(i+1));
    m4 = g*I(i+1) - d*(R(i)+h*m3) + u(i+1)*S(i+1);
    R(i+1) = R(i) + (h/6)*(m1 + 2*m2 + 2*m3 + m4);
end
%Una vez que se obtiene la convergencia deseada, los valores de los vectores finales se almacenan en la matriz 'y'
y(1,:) = t;
y(2,:) = S;
y(3,:) = E;
y(4,:) = I;
y(5,:) = R;

y(6,:) = N;
y(7,:) = u;

y(8,:) = Sn;
y(9,:) = En;
y(10,:) = In;
y(11,:) = Rn;

%Generacion de las graficas

%T = table(S,u);
%disp(t);
%plot(y(1,:),y(2,:)) 
%title("Susceptibles")

%plot(y(1,:),y(6,:)) 
%title("Poblacion Total")

%plot(y(1,:),y(4,:)) 
%title("Infectados")

%plot(y(1,:),y(5,:)) 
%title("Recuperados")

%plot(y(1,:),y(6,:)) 
%title("Todos")

%plot(y(1,:),y(7,:)) 
%title("Control u(t)")

%subplot(3,2,1),plot(y(1,:),y(2,:)),title("Susceptibles");
%subplot(3,2,2),plot(y(1,:),y(3,:)),title("Expuestos");
%subplot(3,2,3),plot(y(1,:),y(4,:)),title("Infectados");
%subplot(3,2,4),plot(y(1,:),y(5,:)),title("Recuperados");
%subplot(3,2,5),plot(y(1,:),y(6,:)),title("Poblacion TOTAL");
%subplot(3,2,6),plot(y(1,:),y(7,:)),title("Control");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold
%%%%%%plot(y(1,:),y(11,:));
plot(y(1,:),y(4,:),y(1,:),y(10,:));
title('Poblacion de Infectados')
%ax = gca;
%ax.FontSize = 18;
legend({'y = control','y = sin control'});%,'Location','northeast')
%hold


%plot(y(1,:),y(2,:),y(1,:),y(8,:));
%title('Poblacion de Susceptible')
%ax = gca;
%ax.FontSize = 18;
%legend({'f = control','f = sin control'},'Location','northeast')

%hold
%plot(y(1,:),y(4,:));
