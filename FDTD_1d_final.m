close all;
clear all;

ic = 80; 
iw = 4; %for a and c
%iw = 30; %forb
yeta = 377; Cf = 0.5;

L = 0.5; N = 500; 
timestep_number = 700; 
dz = 1e-3; dt = 1.67* 1e-12; 
eps0 = 8.854e-12; mu0 = pi*4e-7; c = 1/(eps0.*mu0).^0.5;
zc = ic*dz; w = iw*dz;
x1 = linspace(0,L,N+1); 
x2 = linspace(0,L,N);

cb_dz = ones(N,1)*dt/(dz*eps0);
db_dz = ones(N+1,1)*dt/(dz*mu0);

Ex = zeros(N+1,1);
Hy = zeros(N,1);
Capture_Ex = zeros(N+1,8);
Exact_Ex = zeros(N+1,8);
rotate = 1;
Pulse = zeros(1,N); 

for k = 1:8
    for i = 1:N+1
        z = i * dz; t = (k-1)*100*dt;
        Exact_Ex(i,k) = exp(-((z - zc - c*t)/w)^2); 
    end
    Exact_Ex(1,k) = 0;
    Exact_Ex(N+1,k) = 0;
end
%initial conditions
for k = 1:N
    Pulse(k) = exp(-((k-ic)/iw)^2);
    %Hy(k) = exp(-((k+1-ic-Cf/2)/iw)^2) ./ yeta; %for a and b
    Hy(k) = 0; %for c
end
Ex(1) = 0; Ex(N+1) = 0; 

for iter=0:timestep_number
    
   for k=2:N    %update E
      Ex(k) = Ex(k)-cb_dz(k)*(Hy(k)-Hy(k-1));
   end
    
   for k=1:N   %update H
      Hy(k) = Hy(k)-db_dz(k)*(Ex(k+1)-Ex(k));
   end
  
  
   
  %%%%%%%%%%%%%%%%to make plot capture%%%%%%%%%%%%%%%%%%%%%% 
   if(rem(iter,100) == 0)
       for k = 1 : N-1
        Capture_Ex(k,rotate) = Ex(k);
       end
      rotate = rotate +1;
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%to make the plot move%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %hold off
   %plot(x1,Ex,'b');
   %axis([dz 0.5 -1 1]);
   %grid on;
   %title('Ex');
   %hold on
   %pause(0);
   %iter
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for i = 1 : 8
    subplot(8,1,i)
    plot(x1,Capture_Ex(:,i),'b');
    hold on;
    plot(x1,Exact_Ex(:,i),'r');
    axis([dz 0.5 -1 1]);
end



