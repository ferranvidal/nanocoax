function [mesh,setup,Tp,FE,EDG,HDG,JDG,RDG] = run_example(G,b,lambda)


%% Load mesh and simulation setup

load(['mesh_G=',num2str(G),'.mat'])

f = 300/lambda;
omegapE = 8.45; %eV
gammaE = 0.047; %eV
vf = 1.39e6; %m/s
omegap = omegapE*setup.nondimFactor; % adimensional
gamma = gammaE*setup.nondimFactor; % adimensional
beta2 = b*3/5*(vf/setup.c)^2; % adimensional
omega = f*setup.omega0; % adimensional
setup.tau_n = omegap/sqrt(beta2);
setup.tau_t = f*setup.tau_0;
setup.param = {omega,omegap,gamma,beta2};

%% Materials

if b == 0
    permittivity = [sapphire(lambda),gold(lambda),1,alumina(lambda)];
else
    permittivity = [sapphire(lambda),1,1,alumina(lambda)];
end

k = 1;
for n = 1:mesh.num_mat
    mesh.permittivity(:,mesh.mate == n,:) = permittivity(k);
    k = k + 1;
end


%% Solver
[EDG,HDG,JDG,RDG] = solver_maxwell(mesh,setup);

%% Incident wave
k = setup.k;
p = setup.p;
kxp = cross(k,p);
kX = mesh.dgnodes(:,1,:)*k(1) + ...
    mesh.dgnodes(:,2,:)*k(2) + ...
    mesh.dgnodes(:,3,:)*k(3);

n_incident = sqrt(permittivity(1));
omegaeff = 1i*omega*n_incident;

Ei = zeros(mesh.npv,mesh.nd,mesh.ne);
Hi = zeros(mesh.npv,mesh.nd,mesh.ne);
Ei(:,logical(p),:) = exp(omegaeff*kX);
Hi(:,logical(kxp),:) = n_incident*exp(omegaeff*kX);

%% Transmission
Pt = 0.5*transmission(mesh,EDG,HDG,mesh.power);
P0 = 0.5*transmission(mesh,Ei,Hi,mesh.power0);
Tp = 100*Pt/P0;

%% Field enhancement
[F,V] = field_enhancement(mesh,EDG,mesh.fieldEnh);
FE = F/V;