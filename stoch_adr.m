%THIS CODE SOLVES THE FULL STOCHASTIC ADR SYSTEM

rng(1)
d = 0.01; %dx
h= 5; %dt

Tmax = 10000; %sec
L=1;
n=floor(L/d);

x = linspace(-L, L, n);
y = x;
[X,Y] = meshgrid(x, y);


D = 4.9*10^-6; %diffusione m^2/s ; %diffusione m^2/s
b=1*10^(-9); %era -2
a0=1/1200; %era -0

%DEFINITION OF THE NOISE
m=20; %modes
nu=D/2;
S=zeros(2, n, n, 10, 10);
Ntheta2=0;
for k1 = -m : m
    for  k2 = -m : m
        if(norm([k1, k2])<=floor(0) || norm([k1, k2])>m) 
            theta=0;
        else   
            theta= sqrt(2)/norm([k1, k2])^0;
        end    
            Ntheta2 = Ntheta2 + theta^2;
                if(k1>0 || (k1==0 && k2>0))
                    S(1, :, :, m+1+k1, m+1+k2)= theta*(k2/norm([k1, k2])).*cos(2*pi*(k1*X + k2*Y )); 
                    S(2, :, :, m+1+k1, m+1+k2)= theta*(-k1/norm([k1, k2])).*cos(2*pi*(k1*X + k2*Y )); 
                elseif(k1<0  || (k1==0 && k2<0))
                    S(1, :, :, m+1+k1, m+1+k2)= theta*(k2/norm([k1, k2])).*sin(2*pi*(k1*X + k2*Y )); 
                    S(2, :, :, m+1+k1, m+1+k2)= theta*(-k1/norm([k1, k2])).*sin(2*pi*(k1*X + k2*Y ));
                end
    end
end
eps=(2*sqrt(nu))/sqrt(Ntheta2);
%Ssum(:, :, :)= eps* sum(S, [4,5]);

%CONSTANT WIND AS MULTIPLE OF FISHER SPEED
VF = 2*sqrt(D*a0); %Fisher speed
theta = pi/3;
U = [cos(theta), sin(theta)] * 6 * VF; 


%DISCRETIZING THE LAPLACIAN
A = (D+nu) * spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n) / d^2; %derivata parziale seconda in x
A(1,n)=(D+nu)/d^2; A(n,1)=(D+nu)/d^2;
AD = kron(speye(n), A) + kron(A, speye(n)); %derivata parziale seconda in x + derivata parziale in y
%AA è costruita in modo che agisca su un vettore n^2x1 ottenuto
%concatenando le colonne del mesh

%DISCRETIZING THE GRADIENT
G = spdiags(ones(n, 1) * [-1 0 1], -1:1, n, n) / (2*d); %derivata parziale in x
G(1,n)= -1/(2*d); G(n,1) = 1/(2*d);
AA = U(2) * kron(speye(n), G) + U(1) * kron(G, speye(n)); %U
%grad = [kron(G, speye(n)); kron(speye(n), G)];

%DISCRETIZING THE REACTION
Delta=0.1*a0; %era 10
ran_r = a0 + Delta*(rand(n)-1/2);
        %ran_r= a0 + zeros(n);
        %ran_r(50, 50)= a0 + 1*Delta;
cc=ran_r/b;
C = reshape(ran_r, n^2, 1);
chem = @(F)(ran_r.*F - b*F.^2);  

%Using implicit Euler for advection and diffusion. Esplicit euler
%for the reaction partù

sprintf("CFLs: delta=%f, alpha=%f, rho=%f ", (D+nu)*h/(d^2), norm(U)*h/d, a0*h)
sprintf("Pe(A/D)=%f, DaD(C/D)=%f, DaA(C/A)=%f ", L*(norm(U)/D), a0*2*L^2/D, a0*L/norm(U))
% Risolvo il sistema (I - h/2 * L_D) u_{n+1} = (I + h/2*L_D)(I + h*L_A + sqrt(h)N*L_AR)u_n 
u0=(a0/b)*ones(n);
%u0 = exp(-1/(0.3^2)*((X-0.5).^2 + (Y-0.3).^2));
u=u0;
T = (speye(size(AD)) - h/2 * (AD+AA)) \ (speye(size(AD)) + h/2 * (AD+AA));
gradux = kron(G, speye(n));
graduy = kron(speye(n), G);
phi=zeros(n);
absphi=zeros(n);
figure()
for j = 1 : Tmax/h
    NOISE = zeros(2,n);
    L_Cu=(reshape(u, (n)^2, 1) + h * reshape(chem(u), (n)^2, 1));
        for k1 = 1 : 2*m
            for  k2 = 1 : 2*m
             NOISE = NOISE + normrnd(0,1)*(S(:,:,:, k1, k1)) ;
            end
        end
    NOISE=reshape(NOISE, 2, n^2, 1);
    AN = gradux.* (NOISE(1, : ))' + graduy.* (NOISE(2, : ))' ; 
    unew = T *((speye(n^2)+ eps*sqrt(h)*AN)*L_Cu);
    unew = reshape(unew, n, n);
   

    % DYNAMIC PLOTTING OF THE SOLUTION
    mesh(x, y, unew, 'FaceColor', 'flat');
    view(0,90) %view of the plot from above
    title(sprintf('Time = %f, 0VF', j * h));
    %axis([-1 1 -1 1 a0/b-0.1 a0/b+0.1])
drawnow;
    pause(0.01);
    if (norm(u-unew) < 8*10^-5)
        sprintf("Steady State reached")
        u=unew;
        break
    end
  u=unew;
  phi=phi+fft2(b*(u-u0)/a0);
  absphi=absphi + abs(fftshift(phi));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTING AND PLOTTING OF THE CORRELATION FUNCTIONS

phi=phi/j;
absphi=absphi/j;
figure()
imagesc(log(absphi))
colorbar
title(sprintf('E[|phi(q)|^2], 3VF, Noise'));

%correlation function
corr=zeros(n);
u_=normalize(u);
for i=1:n
    for j=1:n
        corr(i, j)=mean( u_.*circshift(u_, [i,j]),'all');
    end
end
 mesh(x, y, circshift(corr, [50,50]), 'FaceColor', 'flat');
    view(0,90) %view of the plot from above

