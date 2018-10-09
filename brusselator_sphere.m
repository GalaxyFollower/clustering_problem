%% Reaction-diffusion on a sphere

%cpf = @cpSphere;
%paramf = @paramSphere;
%cpf = @cpEllipsoid;
%paramf = @paramEllipsoid;

R = 1.0;          % cap diameter

cpf = @(x,y,z) cpSphere(x, y, z, R);
paramf = @(N) paramSphere(N, R);


dx = 0.05;      % grid size default = 0.05

% make vectors of x, y, z positions of the grid
x1d = ((-R-5*dx):dx:(R+5*dx))';
y1d = x1d;
z1d = x1d;

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

[cpx, cpy, cpz, dist] = cpf(xx,yy,zz);
%[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz,R);

%% Banding: do calculation in a narrow band around the surface
dim = 3;  % dimension
p = 3;    % interpolation order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
x = xx(band); y = yy(band); z = zz(band);

%% discrete operators
disp('building laplacian and interp matrices');
L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
E1 = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, 1, band);
E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
I = speye(size(E));

%% plotting grid
[xp,yp,zp] = paramf(64);
Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);

figure(1); clf;
figure(2); clf;

% u_t = f(u,g) + nuu*Lap u
% v_t = g(u,g) + nuv*Lap u

% parameters from [Ward Ruuth, todo, 2014]
epsilon = 0.075;
%epsilon = 0.1;
EE = 4;
ff = 0.8;
tau = 1/ff^2
nonlinu = @(u,v) epsilon^2*EE - u + ff*(u.^2.*v);
nonlinv = @(u,v) (1/(tau*epsilon^2))*(u - u.^2.*v);

Du = epsilon^2;
Dv = 0.2/tau;

ue = epsilon^2*EE/(1-ff)
ve = (1-ff)/(EE*epsilon^2)



%% initial conditions - perturbation from steady state
%pert = 0.5*exp(-(10*(z-.1)).^2) + 0.5*sin(5*pi*x) + 0.5*rand(size(x));
%u0 = ue*(1 + 0.1*pert);
%v0 = ve*(1 + 0.1*pert);

% perturbation of 2% from special ue, ve values
u0 = ue*(1 + 0.02*rand(size(x)));
v0 = ve*(1 + 0.02*rand(size(x)));
u = u0;  v = v0;


%% time-stepping
Tf = 70;
implicit = true;
if (implicit)
  dt = 1 * (1/(tau*Dv)) * dx
else
  dt1 = 1/8 * (1/max(Du,Dv)) * dx^2
  dt = 1/(8*tau*Dv)*dx^2
end
numtimesteps = ceil(Tf/dt)
% adjust for integer number of steps
dt = Tf / numtimesteps


figure(1);
sphplot = Eplot*u;
sphplot = reshape(sphplot, size(xp));
Hplot = surf(xp, yp, zp, sphplot);
title('initial u')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
%axis off;
shading interp
%camlight left
colorbar

figure(2);
sphplot = Eplot*v;
sphplot = reshape(sphplot, size(xp));
Hplot2 = surf(xp, yp, zp, sphplot);
title('initial v')
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
view(-10, 60)
%axis off;
shading interp
%camlight left
colorbar


%% Method-of-lines approach
% See [vonGlehn/Macdonald/Maerz 2013]
lambda = 6/(dx^2);
%lambdav = 6*Dv/(dx^2);
LBu = Du*((E1*L) - lambda*(I-E));
LBv = Dv*((E1*L) - lambda*(I-E));

Au = I - dt*LBu;
Av = I - dt*LBv;

plotskip = min(ceil(0.5/dt), 50)

tic; [pcLu, pcUu] = ilu (Au); toc
tic; [pcLv, pcUv] = ilu (Av); toc

e1 = tic();
for kt = 1:numtimesteps
  if (~ implicit)
    %% MOL: explicit Euler timestepping
    %unew = u + dt*( E*nonlinu(u,v) + LBu*u );
    %vnew = v + dt*( E*nonlinv(u,v) + LBv*v );
    %unew = u + dt*( nonlinu(u,v) + LBu*u );
    %vnew = v + dt*( nonlinv(u,v) + LBv*v );
    %% RM iteration
    unew = u + dt*( nonlinu(u,v) + Du*(L*u));
    vnew = v + dt*( nonlinv(u,v) + Dv*(L*v));
    unew = E*unew;
    vnew = E*vnew;
  else
    %% MOL: implicit Euler
    tic
    %unew = Au \ (u + dt*(nonlinu(u,v)));
    %vnew = Av \ (v + dt*(nonlinv(u,v)));
    unew = gmres(Au, (u + dt*(nonlinu(u,v))), 12, 1e-8, [], pcLu, pcUu);
    toc
    tic
    vnew = gmres(Av, (v + dt*(nonlinv(u,v))), 12, 1e-8, [], pcLv, pcUv);
    toc
  end
  u = unew;
  v = vnew;
  %% MOL: without precomputing matrices
  %rhsu = Du*(L*u) + nonlinu(u,v);
  %rhsv = Dv*(L*v) + nonlinv(u,v);
  %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
  %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );

  %% Ruuth-Merriman
  %rhsu = nuu*(L*u) + f(u,v);
  %rhsv = nuv*(L*v) + g(u,v);
  %unew = u + dt*rhsu;
  %vnew = v + dt*rhsv;
  %u = E*unew;
  %v = E*vnew;

  t = kt*dt;

  if ( (mod(kt,plotskip)==0) || (kt<=20) || (kt==numtimesteps) )
    disp([kt t min(u) max(u) min(v) max(v)]);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 1);
    set(Hplot, 'CData', sphplot);
    title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
    sphplot = Eplot*v;
    sphplot = reshape(sphplot, size(xp));
    set(0, 'CurrentFigure', 2);
    set(Hplot2, 'CData', sphplot);
    drawnow;
  end
end
toc(e1)
