eps0 = 10.^-9/36/pi; mu0 = 4*pi*10.^-7;
c = 3*10.^8; Cn = 0.99;
dx = 0.001; dy =0.001; dz = 0.001;
dt = Cn/(c * ((1/dz).^2+(1/dy).^2+(1/dz).^2).^.5); 
dt_mu = dt/mu0; dt_eps = dt/eps0;

duration = 8192;
nx = 30; ny = 40; nz = 10;

Ex = zeros(nx,ny+1,nz+1); Hx = zeros(nx+1,ny,nz);
Ey = zeros(nx+1,ny,nz+1); Hy = zeros(nx,ny+1,nz);
Ez = zeros(nx+1,ny+1,nz); Hz = zeros(nx,ny,nz+1);
Ez_save = zeros(duration,1);

Pulse = zeros(duration,1);
T = 2.^.5*(log(2)).^.5 /(pi*15*10^9);
for n = 1:duration
    Pulse(n) = exp(-((n*dt - 3*T)/T).^2);
end

for iter = 1:duration
 %%%%%%%%%%%Ecell update%%%%%%%%%%%%%%%%%%%%%%
 Ex(2:nx,2:ny,2:nz) = Ex(2:nx,2:ny,2:nz)...
               +(dt_eps/dy).*(Hz(2:nx,2:ny,2:nz)-Hz(2:nx,1:ny-1,2:nz))...
               -(dt_eps/dz).*(Hy(2:nx,2:ny,2:nz)-Hy(2:nx,2:ny,1:nz-1));

Ey(2:nx,2:ny,2:nz) = Ey(2:nx,2:ny,2:nz)...
               +(dt_eps/dz)*(Hx(2:nx,2:ny,2:nz)-Hx(2:nx,2:ny,1:nz-1))...
               -(dt_eps/dx)*(Hz(2:nx,2:ny,2:nz)-Hz(1:nx-1,2:ny,2:nz));                              

Ez(2:nx,2:ny,2:nz) = Ez(2:nx,2:ny,2:nz)...
               +(dt_eps/dx)*(Hy(2:nx,2:ny,2:nz)-Hy(1:nx-1,2:ny,2:nz))...
               -(dt_eps/dy)*(Hx(2:nx,2:ny,2:nz)-Hx(2:nx,1:ny-1,2:nz));
   
   %%%%%%%%%%%%%source 
   Ez(3,3,3) =  Ez(3,3,3) + Pulse(iter);
  
   
  %%%%%%%%%Hcell update
 Hx(1:nx+1,1:ny,1:nz) =Hx(1:nx+1,1:ny,1:nz)...
                 +(dt_mu/dz)*(Ey(1:nx+1,1:ny,2:nz+1)-Ey(1:nx+1,1:ny,1:nz))...
                 -(dt_mu/dy)*(Ez(1:nx+1,2:ny+1,1:nz)-Ez(1:nx+1,1:ny,1:nz));

Hy(1:nx,1:ny+1,1:nz) = Hy(1:nx,1:ny+1,1:nz)...
                 +(dt_mu/dx)*(Ez(2:nx+1,1:ny+1,1:nz)-Ez(1:nx,1:ny+1,1:nz))...
                 -(dt_mu/dz)*(Ex(1:nx,1:ny+1,2:nz+1)-Ex(1:nx,1:ny+1,1:nz));             
                          
Hz(1:nx,1:ny,1:nz+1) = Hz(1:nx,1:ny,1:nz+1)...
                 +(dt_mu/dy)*(Ex(1:nx,2:ny+1,1:nz+1)-Ex(1:nx,1:ny,1:nz+1))...
                 -(dt_mu/dx)*(Ey(2:nx+1,1:ny,1:nz+1)-Ey(1:nx,1:ny,1:nz+1));  
   
  Ez_save(iter) = Ez(12,3,10);
   
end
subplot(2,1,1);
plot(1:duration,Ez_save);

subplot(2,1,2);
[fax, fdata] = dtft(Ez_save, dt, 5*1e9, 30*1e9, 1e7);
mag_Ez = abs(fdata);
plot(fax, mag_Ez);
xticks([5*1e9:5*1e9:30*1e9]);
xticklabels({'5', '10', '15', '20', '25', '30'});

