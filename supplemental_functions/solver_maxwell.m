function [EDG,HDG,JDG,RDG]= solver_maxwell(mesh,setup)


%% Setup of parameters, elemental products of shape functions, Jacobians and normals

master = mesh.master;
nd  = mesh.nd;
ne  = mesh.ne;
npv = master.npv;
ngv = master.ngv;
ngf = master.ngf;
npf = master.npf;
nfe = size(master.perm,2);
nf = mesh.nf;
nch = nd-1;
ncu = nd;
ncv = nd;
ncj = nd;
ncq = 1;
ncqh = 1;
ncht = nch + ncqh;
mesh.metal = reshape(mesh.metal,1,[]);


omega = setup.param{1};
omegap = setup.param{2};
gamma = setup.param{3};
beta2 = setup.param{4};
if beta2 == 0
    ishydro = false;
else
    ishydro = true;
end
alpha = omega*(omega+1i*gamma);
delta = 1i*omega*omegap^2;
theta = 1i*omega;
tau_t = setup.tau_t;
tau_n = setup.tau_n;
bcm_t  = setup.bcm_t;
bcm_n = setup.bcm_n;
bf   = mesh.bf;


mu = mesh.permeability;
ep = mesh.permittivity;


shapvt = master.shapvt(:,:,1);
dshapvt = reshape(permute(master.shapvt(:,:,2:end),[1 3 2]),[ngv*nd npv]);
shapfc = squeeze(master.shapfc(:,:,1));
shapft  = master.shapfc(:,:,1)';
dshapft = reshape(permute(master.shapft(:,:,2:end),[1 3 2]),[ngf*(nd-1) npf]);

perm = master.perm;
gwfc = master.gwfc;


shapvgdotshapvl = reshape(bsxfun(@times,reshape(master.shapvg(:,:,1),[npv 1 ngv]),reshape(master.shapvl(:,:,1),[1 npv ngv])),[npv*npv ngv]);
shapfgdotshapfc(:,:,1)  = reshape(bsxfun(@times,reshape(master.shapfg(:,:,1),[npf 1 ngf]),reshape(master.shapfc(:,:,1),[1 npf ngf])),[npf*npf ngf]);
for d = 1:nd
    shapvgdotdshapvl(:,:,d) = reshape(bsxfun(@times,reshape(master.shapvg(:,:,d+1),[npv 1 ngv]),reshape(master.shapvl(:,:,1),[1 npv ngv])),[npv*npv ngv]);   
end

pdg = reshape(mesh.dgnodes(:,1:nd,:),[npv nd*ne]);
Jg = permute(reshape(dshapvt*pdg,[ngv nd nd ne]),[1 4 2 3]);

detJ = Jg(:,:,1,1).*Jg(:,:,2,2).*Jg(:,:,3,3) - Jg(:,:,1,1).*Jg(:,:,3,2).*Jg(:,:,2,3)+ ...
    Jg(:,:,2,1).*Jg(:,:,3,2).*Jg(:,:,1,3) - Jg(:,:,2,1).*Jg(:,:,1,2).*Jg(:,:,3,3)+ ...
    Jg(:,:,3,1).*Jg(:,:,1,2).*Jg(:,:,2,3) - Jg(:,:,3,1).*Jg(:,:,2,2).*Jg(:,:,1,3);
Jinv11 = (Jg(:,:,2,2).*Jg(:,:,3,3) - Jg(:,:,2,3).*Jg(:,:,3,2));
Jinv21 = (Jg(:,:,2,3).*Jg(:,:,3,1) - Jg(:,:,3,3).*Jg(:,:,2,1));
Jinv31 = (Jg(:,:,2,1).*Jg(:,:,3,2) - Jg(:,:,2,2).*Jg(:,:,3,1));
Jinv12 = (Jg(:,:,3,2).*Jg(:,:,1,3) - Jg(:,:,3,3).*Jg(:,:,1,2));
Jinv22 = (Jg(:,:,1,1).*Jg(:,:,3,3) - Jg(:,:,1,3).*Jg(:,:,3,1));
Jinv32 = (Jg(:,:,1,2).*Jg(:,:,3,1) - Jg(:,:,3,2).*Jg(:,:,1,1));
Jinv13 = (Jg(:,:,1,2).*Jg(:,:,2,3) - Jg(:,:,2,2).*Jg(:,:,1,3));
Jinv23 = (Jg(:,:,1,3).*Jg(:,:,2,1) - Jg(:,:,1,1).*Jg(:,:,2,3));
Jinv33 = (Jg(:,:,1,1).*Jg(:,:,2,2) - Jg(:,:,1,2).*Jg(:,:,2,1));

pfg = permute(reshape(shapft*reshape(pdg(perm,:),[npf nfe*ne*nd]),[ngf nfe nd ne]),[1 3 4 2]);
dpfg = permute(reshape(dshapft*reshape(pdg(perm,:),[npf nfe*ne*nd]),[ngf nd-1 nfe nd ne]),[1 5 3 4 2]);

Nx = dpfg(:,:,:,2,1).*dpfg(:,:,:,3,2) - dpfg(:,:,:,3,1).*dpfg(:,:,:,2,2);
Ny = dpfg(:,:,:,3,1).*dpfg(:,:,:,1,2) - dpfg(:,:,:,1,1).*dpfg(:,:,:,3,2);
Nz = dpfg(:,:,:,1,1).*dpfg(:,:,:,2,2) - dpfg(:,:,:,2,1).*dpfg(:,:,:,1,2);
detJf = sqrt(Nx.^2 + Ny.^2 + Nz.^2);
Nx = Nx./detJf; Ny = Ny./detJf; Nz = Nz./detJf;

if ishydro
    eleL = find(mesh.metal==0);
    eleNL = find(mesh.metal==1);
    neL  = numel(eleL);
    neNL  = numel(eleNL);
else
    eleL = 1:ne;
    neL = ne;
    eleNL = [];
    neNL = 0;
end

warning('off','MATLAB:nearlySingularMatrix')

%% DIELECTRIC, solving Maxwell's equations
tic
WL  = zeros(nfe*npf*nch,nfe*npf*nch,neL);
ZL  = zeros((ncu+ncv)*npv,nfe*npf*nch,neL);
fprintf('Local elements...')
if neL > 0
    
    ngrsiz = 256;
    ngr    = ceil(neL/ngrsiz);
    ngrne  = round(neL/ngr);
    nk = 1:ngrne:neL;
    nb = [nk(1:end); [nk(2:end)-1,neL]];
    blocksL = nb;
    % Run assembly in chunks of elements
    
    for b = 1:size(blocksL,2)
        e1 = blocksL(1,b);
        e2 = blocksL(2,b);
        eleLb = eleL(e1:e2);
        neLb = (e2-e1+1);
        
        mug = shapvt*mu(:,eleLb);
        epxg = shapvt*ep(:,eleLb,1);
        epyg = shapvt*ep(:,eleLb,2);
        epzg = shapvt*ep(:,eleLb,3);
        
        % Volume matrices
        AAmu = zeros(npv,npv,neLb,ncv);
        AAep = zeros(npv,npv,neLb,ncu);
        BB = zeros(npv,npv,neLb,nd);
       
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*mug);
        AAmu(:,:,:,1) = reshape(tm,[npv npv neLb]);
        
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*mug);
        AAmu(:,:,:,2) = reshape(tm,[npv npv neLb]);
        
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*mug);
        AAmu(:,:,:,3) = reshape(tm,[npv npv neLb]);
        
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*epxg);
        AAep(:,:,:,1) = reshape(tm,[npv npv neLb]);
        
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*epyg);
        AAep(:,:,:,2) = reshape(tm,[npv npv neLb]);
        
        tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleLb).*epzg);
        AAep(:,:,:,3) = reshape(tm,[npv npv neLb]);
        
        tmx = reshape(shapvgdotdshapvl(:,:,1)*Jinv11(:,eleLb) +...
            shapvgdotdshapvl(:,:,2)*Jinv12(:,eleLb) + shapvgdotdshapvl(:,:,3)*Jinv13(:,eleLb),[npv npv neLb]);
        BB(:,:,:,1) = tmx;
        
        tmy = reshape(shapvgdotdshapvl(:,:,1)*Jinv21(:,eleLb) +...
            shapvgdotdshapvl(:,:,2)*Jinv22(:,eleLb) + shapvgdotdshapvl(:,:,3)*Jinv23(:,eleLb),[npv npv neLb]);
        
        BB(:,:,:,2) = tmy;
        
        tmz = reshape(shapvgdotdshapvl(:,:,1)*Jinv31(:,eleLb) +...
            shapvgdotdshapvl(:,:,2)*Jinv32(:,eleLb) + shapvgdotdshapvl(:,:,3)*Jinv33(:,eleLb),[npv npv neLb]);
        
        BB(:,:,:,3) = tmz;

        % Face matrices
        DD = zeros(npv,npv,neLb,ncu*(ncu+1)/2);
        EE = zeros(npv,npf*nfe,neLb,ncu*nch);
        CC = zeros(npv,npf*nfe,neLb,ncv*nch);
        LL = zeros(npf*nfe,npv,neLb,ncu*nch);
        RR = zeros(npf*nfe,npv,neLb,ncv*nch);
        MM = zeros(npf*nfe,npf*nfe,neLb,nch*(nch+1)/2);
        for jf=1:nfe
            I = perm(:,jf);
            J = ((jf-1)*npf+1):jf*npf;
            
            nx = Nx(:,eleLb,jf);
            ny = Ny(:,eleLb,jf);
            nz = Nz(:,eleLb,jf);
            [t1x,t1y,t1z,t2x,t2y,t2z] = tangentvectors(nx,ny,nz);
            detjf = detJf(:,eleLb,jf);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1x,[ngf neLb]);
            EE(I,J,:,1) = EE(I,J,:,1) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1y,[ngf neLb]);
            EE(I,J,:,2) = EE(I,J,:,2) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1z,[ngf neLb]);
            EE(I,J,:,3) = EE(I,J,:,3) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2x,[ngf neLb]);
            EE(I,J,:,4) = EE(I,J,:,4) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2y,[ngf neLb]);
            EE(I,J,:,5) = EE(I,J,:,5) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2z,[ngf neLb]);
            EE(I,J,:,6) = EE(I,J,:,6) + reshape(tm,[npf npf neLb]);
            

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*t1z - nz.*t1y),[ngf neLb]);
            CC(I,J,:,1) = CC(I,J,:,1) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nz.*t1x - nx.*t1z),[ngf neLb]);
            CC(I,J,:,2) = CC(I,J,:,2) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*t1y - ny.*t1x),[ngf neLb]);
            CC(I,J,:,3) = CC(I,J,:,3) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*t2z - nz.*t2y),[ngf neLb]);
            CC(I,J,:,4) = CC(I,J,:,4) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nz.*t2x - nx.*t2z),[ngf neLb]);
            CC(I,J,:,5) = CC(I,J,:,5) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*t2y - ny.*t2x),[ngf neLb]);
            CC(I,J,:,6) = CC(I,J,:,6) + reshape(tm,[npf npf neLb]);
            

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*ny+nz.*nz),[ngf neLb]);
            DD(I,I,:,1) = DD(I,I,:,1) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*nx+nz.*nz),[ngf neLb]);
            DD(I,I,:,2) = DD(I,I,:,2) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*ny+nx.*nx),[ngf neLb]);
            DD(I,I,:,3) = DD(I,I,:,3) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*ny),[ngf neLb]);
            DD(I,I,:,4) = DD(I,I,:,4) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*nz),[ngf neLb]);
            DD(I,I,:,5) = DD(I,I,:,5) + reshape(tm,[npf npf neLb]);

            tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*nz),[ngf neLb]);
            DD(I,I,:,6) = DD(I,I,:,6) + reshape(tm,[npf npf neLb]);
            
            
            ep_ = ep(:,eleLb,:);
            bfele = bf(jf,eleLb);
            idbou = find(bfele<0);
            id1 = bcm_t(-bfele(idbou))==1; 
            id3 = bcm_t(-bfele(idbou))==3; 
            ef = shapft*ep_(I,idbou(id3),1);
            
            
            const = ones(ngf,neLb);
            const(:,idbou(id1)) = 1/tau_t;
            const(:,idbou(id3)) = 1 - 1i*omega*sqrt(ef)/tau_t;
            cdetjf = const.*detjf;
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t1x.*t1x+t1y.*t1y+t1z.*t1z),[ngf neLb]);
            MM(J,J,:,1) = MM(J,J,:,1) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t1x.*t2x+t1y.*t2y+t1z.*t2z),[ngf neLb]);
            MM(J,J,:,2) = MM(J,J,:,2) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t2x.*t2x+t2y.*t2y+t2z.*t2z),[ngf neLb]);
            MM(J,J,:,3) = MM(J,J,:,3) + reshape(tm,[npf npf neLb]);
            
            
            const = ones(ngf,neLb);
            const(:,idbou(id1)) = 0;
            cdetjf = const.*detjf;
            
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-ny.*t1z + nz.*t1y),[ngf neLb]);
            RR(J,I,:,1) = RR(J,I,:,1) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nz.*t1x + nx.*t1z),[ngf neLb]);
            RR(J,I,:,2) = RR(J,I,:,2) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nx.*t1y + ny.*t1x),[ngf neLb]);
            RR(J,I,:,3) = RR(J,I,:,3) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-ny.*t2z + nz.*t2y),[ngf neLb]);
            RR(J,I,:,4) = RR(J,I,:,4) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nz.*t2x + nx.*t2z),[ngf neLb]);
            RR(J,I,:,5) = RR(J,I,:,5) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nx.*t2y + ny.*t2x),[ngf neLb]);
            RR(J,I,:,6) = RR(J,I,:,6) + reshape(tm,[npf npf neLb]);
            
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1x,[ngf neLb]);
            LL(J,I,:,1) = LL(J,I,:,1) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1y,[ngf neLb]);
            LL(J,I,:,2) = LL(J,I,:,2) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1z,[ngf neLb]);
            LL(J,I,:,3) = LL(J,I,:,3) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2x,[ngf neLb]);
            LL(J,I,:,4) = LL(J,I,:,4) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2y,[ngf neLb]);
            LL(J,I,:,5) = LL(J,I,:,5) + reshape(tm,[npf npf neLb]);
            
            tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2z,[ngf neLb]);
            LL(J,I,:,6) = LL(J,I,:,6) + reshape(tm,[npf npf neLb]);

        end
        
        wL  = zeros(nfe*npf*nch,nfe*npf*nch,neLb);
        zL  = zeros((ncu+ncv)*npv,nfe*npf*nch,neLb);
        % Compute global quantities to assemble element by element
        for ie=1:neLb
            
            zv = zeros(npv);
            
            Y1 = ([AAmu(:,:,ie,1) zv zv; ...
                zv AAmu(:,:,ie,2) zv; ...
                zv zv AAmu(:,:,ie,3)]);
            
            B_K = [zv BB(:,:,ie,3) -BB(:,:,ie,2); ...
                -BB(:,:,ie,3) zv BB(:,:,ie,1); ...
                BB(:,:,ie,2) -BB(:,:,ie,1) zv];
            
            C_K = [CC(:,:,ie,1) CC(:,:,ie,4); ...
                CC(:,:,ie,2) CC(:,:,ie,5); ...
                CC(:,:,ie,3) CC(:,:,ie,6)];
            
            D_K = [DD(:,:,ie,1) -DD(:,:,ie,4) -DD(:,:,ie,5); ...
                -DD(:,:,ie,4)  DD(:,:,ie,2) -DD(:,:,ie,6); ...
                -DD(:,:,ie,5) -DD(:,:,ie,6)  DD(:,:,ie,3)];
            
            D_K = tau_t*D_K - omega^2*[AAep(:,:,ie,1) zv zv; ...
                zv AAep(:,:,ie,2) zv; ...
                zv zv AAep(:,:,ie,3)];
            
            E_K = tau_t*[EE(:,:,ie,1) EE(:,:,ie,4); ...
                EE(:,:,ie,2) EE(:,:,ie,5); ...
                EE(:,:,ie,3) EE(:,:,ie,6)];
            
            
            M_K = tau_t*[MM(:,:,ie,1) MM(:,:,ie,2); ...
                MM(:,:,ie,2) MM(:,:,ie,3)];
            
            L_K = tau_t*[LL(:,:,ie,1) LL(:,:,ie,2) LL(:,:,ie,3); ...
                LL(:,:,ie,4) LL(:,:,ie,5) LL(:,:,ie,6)];
            
            R_K = [RR(:,:,ie,1) RR(:,:,ie,2) RR(:,:,ie,3); ...
                RR(:,:,ie,4) RR(:,:,ie,5) RR(:,:,ie,6)];
   
            Y2 = -B_K;
            Y3 = B_K.';
            Y4 = D_K;
            Y5 = -C_K;
            Y6 = -E_K;
            Y7 = -R_K;
            Y8 = -L_K;
            Y9 = M_K;
            
            YY = [Y1 Y2;Y3 Y4];
            tm = -YY\[Y5;Y6];
            zL(:,:,ie) = tm;
            wL(:,:,ie) = Y9 + [Y7 Y8]*tm;
        end
        ZL(:,:,e1:e2) = zL;
        WL(:,:,e1:e2) = wL;
        clear AAmu AAep BB CC DD EE LL RR MM zL wL
    end
    
end

%%  METAL, solving Maxwell + hydrodynamic model
fprintf('nonlocal elements...')

if neNL > 0
    ngrsiz = 128;
    ngr = ceil(neNL/ngrsiz);
    ngrne  = round(neNL/ngr);
    nk = 1:ngrne:neNL;
    nb = [nk(1:end); [nk(2:end)-1,neNL]];
    blocksNL = nb;
else
    blocksNL = [];
end
    % Run assembly in chunks of elements

WNL  = zeros(nfe*npf*ncht,nfe*npf*ncht,neNL);
ZNL  = zeros((ncu+ncv+ncj+ncq)*npv,nfe*npf*ncht,neNL);

for b = 1:size(blocksNL,2)
    e1 = blocksNL(1,b);
    e2 = blocksNL(2,b);
    eleNLb = eleNL(e1:e2);
    neNLb = (e2-e1+1);
    
    
    mug = shapvt*mu(:,eleNLb);
    epxg = shapvt*ep(:,eleNLb,1);
    epyg = shapvt*ep(:,eleNLb,2);
    epzg = shapvt*ep(:,eleNLb,3);
    
    % Volume matrices
    AAmu = zeros(npv,npv,neNLb,ncv);
    AAep = zeros(npv,npv,neNLb,ncu);
    BB = zeros(npv,npv,neNLb,nd);
    
    AA = reshape(shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb)),[npv npv neNLb]);

    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*mug);
    AAmu(:,:,:,1) = reshape(tm,[npv npv neNLb]);
    
    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*mug);
    AAmu(:,:,:,2) = reshape(tm,[npv npv neNLb]);
    
    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*mug);
    AAmu(:,:,:,3) = reshape(tm,[npv npv neNLb]);
    
    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*epxg);
    AAep(:,:,:,1) = reshape(tm,[npv npv neNLb]);
    
    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*epyg);
    AAep(:,:,:,2) = reshape(tm,[npv npv neNLb]);
    
    tm = shapvgdotshapvl(:,:,1)*(detJ(:,eleNLb).*epzg);
    AAep(:,:,:,3) = reshape(tm,[npv npv neNLb]);
    
    
    tmx = reshape(shapvgdotdshapvl(:,:,1)*Jinv11(:,eleNLb) +...
        shapvgdotdshapvl(:,:,2)*Jinv12(:,eleNLb) + shapvgdotdshapvl(:,:,3)*Jinv13(:,eleNLb),[npv npv neNLb]);
    BB(:,:,:,1) = tmx;
    
    tmy = reshape(shapvgdotdshapvl(:,:,1)*Jinv21(:,eleNLb) +...
        shapvgdotdshapvl(:,:,2)*Jinv22(:,eleNLb) + shapvgdotdshapvl(:,:,3)*Jinv23(:,eleNLb),[npv npv neNLb]);
    
    BB(:,:,:,2) = tmy;
    
    tmz = reshape(shapvgdotdshapvl(:,:,1)*Jinv31(:,eleNLb) +...
        shapvgdotdshapvl(:,:,2)*Jinv32(:,eleNLb) + shapvgdotdshapvl(:,:,3)*Jinv33(:,eleNLb),[npv npv neNLb]);
    
    BB(:,:,:,3) = tmz;

    
    % Face matrices
    DD = zeros(npv,npv,neNLb,ncu*(ncu+1)/2);
    HH = zeros(npv,npv,neNLb);
    NN = zeros(npv,npf*nfe,neNLb);
    MM = zeros(npf*nfe,npf*nfe,neNLb,nch*(nch+1)/2);
    RR = zeros(npf*nfe,npv,neNLb,ncv*nch);
    LL = zeros(npf*nfe,npv,neNLb,ncu*nch);
    EE = zeros(npv,npf*nfe,neNLb,ncu*nch);
    CC = zeros(npv,npf*nfe,neNLb,ncv*nch);
    GG = zeros(npv,npf*nfe,neNLb,ncj);
    TT = zeros(npf*nfe,npf*nfe,neNLb);
    
    for jf=1:nfe
        I = perm(:,jf);
        J = ((jf-1)*npf+1):jf*npf;
        
        nx = Nx(:,eleNLb,jf);
        ny = Ny(:,eleNLb,jf);
            nz = Nz(:,eleNLb,jf);
            [t1x,t1y,t1z,t2x,t2y,t2z] = tangentvectors(nx,ny,nz);
        detjf = detJf(:,eleNLb,jf);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1x,[ngf neNLb]);
        EE(I,J,:,1) = EE(I,J,:,1) + reshape(tm,[npf npf neNLb]);

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1y,[ngf neNLb]);
        EE(I,J,:,2) = EE(I,J,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*nx,[ngf neNLb]);
        GG(I,J,:,1) = GG(I,J,:,1) + reshape(tm,[npf npf neNLb]);
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*ny,[ngf neNLb]);
        GG(I,J,:,2) = GG(I,J,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t1z,[ngf neNLb]);
        EE(I,J,:,3) = EE(I,J,:,3) + reshape(tm,[npf npf neNLb]);

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2x,[ngf neNLb]);
        EE(I,J,:,4) = EE(I,J,:,4) + reshape(tm,[npf npf neNLb]);

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2y,[ngf neNLb]);
        EE(I,J,:,5) = EE(I,J,:,5) + reshape(tm,[npf npf neNLb]);

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*t2z,[ngf neNLb]);
        EE(I,J,:,6) = EE(I,J,:,6) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*nz,[ngf neNLb]);
        GG(I,J,:,3) = GG(I,J,:,3) + reshape(tm,[npf npf neNLb]);
        
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*t1z - nz.*t1y),[ngf neNLb]);
        CC(I,J,:,1) = CC(I,J,:,1) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nz.*t1x - nx.*t1z),[ngf neNLb]);
        CC(I,J,:,2) = CC(I,J,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*t1y - ny.*t1x),[ngf neNLb]);
        CC(I,J,:,3) = CC(I,J,:,3) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*t2z - nz.*t2y),[ngf neNLb]);
        CC(I,J,:,4) = CC(I,J,:,4) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nz.*t2x - nx.*t2z),[ngf neNLb]);
        CC(I,J,:,5) = CC(I,J,:,5) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*t2y - ny.*t2x),[ngf neNLb]);
        CC(I,J,:,6) = CC(I,J,:,6) + reshape(tm,[npf npf neNLb]);
        
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*ny+nz.*nz),[ngf neNLb]);
        DD(I,I,:,1) = DD(I,I,:,1) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*nx+nz.*nz),[ngf neNLb]);
        DD(I,I,:,2) = DD(I,I,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*ny+nx.*nx),[ngf neNLb]);
        DD(I,I,:,3) = DD(I,I,:,3) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*ny),[ngf neNLb]);
        DD(I,I,:,4) = DD(I,I,:,4) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(nx.*nz),[ngf neNLb]);
        DD(I,I,:,5) = DD(I,I,:,5) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detjf.*(ny.*nz),[ngf neNLb]);
        DD(I,I,:,6) = DD(I,I,:,6) + reshape(tm,[npf npf neNLb]);
        
        
        bfeleNL = bf(jf,eleNLb);
        idbou = find(bfeleNL<0);
        

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf,[ngf neNLb]);
        HH(I,I,:) = HH(I,I,:) + reshape(tm,[npf npf neNLb]);
        NN(I,J,:) = NN(I,J,:) + reshape(tm,[npf npf neNLb]);
        

        tm = shapfgdotshapfc(:,:,1)*reshape(detjf,[ngf neNLb]);
        TT(J,J,:) = TT(J,J,:) + reshape(tm,[npf npf neNLb]);        
        
        id1 = bcm_t(-bfeleNL(idbou))==1; 
        
        const = ones(ngf,neNLb);
        const(:,idbou(id1)) = 1/tau_t;
        cdetjf = const.*detjf;
        


        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t1x.*t1x+t1y.*t1y+t1z.*t1z),[ngf neNLb]);
        MM(J,J,:,1) = MM(J,J,:,1) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t1x.*t2x+t1y.*t2y+t1z.*t2z),[ngf neNLb]);
        MM(J,J,:,2) = MM(J,J,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(t2x.*t2x+t2y.*t2y+t2z.*t2z),[ngf neNLb]);
        MM(J,J,:,3) = MM(J,J,:,3) + reshape(tm,[npf npf neNLb]);
        
        
        const = ones(ngf,neNLb);
        const(:,idbou(id1)) = 0;
        cdetjf = const.*detjf;
        

        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-ny.*t1z + nz.*t1y),[ngf neNLb]);
        RR(J,I,:,1) = RR(J,I,:,1) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nz.*t1x + nx.*t1z),[ngf neNLb]);
        RR(J,I,:,2) = RR(J,I,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nx.*t1y + ny.*t1x),[ngf neNLb]);
        RR(J,I,:,3) = RR(J,I,:,3) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-ny.*t2z + nz.*t2y),[ngf neNLb]);
        RR(J,I,:,4) = RR(J,I,:,4) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nz.*t2x + nx.*t2z),[ngf neNLb]);
        RR(J,I,:,5) = RR(J,I,:,5) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*(-nx.*t2y + ny.*t2x),[ngf neNLb]);
        RR(J,I,:,6) = RR(J,I,:,6) + reshape(tm,[npf npf neNLb]);
        
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1x,[ngf neNLb]);
        LL(J,I,:,1) = LL(J,I,:,1) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1y,[ngf neNLb]);
        LL(J,I,:,2) = LL(J,I,:,2) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t1z,[ngf neNLb]);
        LL(J,I,:,3) = LL(J,I,:,3) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2x,[ngf neNLb]);
        LL(J,I,:,4) = LL(J,I,:,4) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2y,[ngf neNLb]);
        LL(J,I,:,5) = LL(J,I,:,5) + reshape(tm,[npf npf neNLb]);
        
        tm = shapfgdotshapfc(:,:,1)*reshape(cdetjf.*t2z,[ngf neNLb]);
        LL(J,I,:,6) = LL(J,I,:,6) + reshape(tm,[npf npf neNLb]);
        
        
    end
    
    wNL  = zeros(nfe*npf*ncht,nfe*npf*ncht,neNLb);
    zNL  = zeros((ncu+ncv+ncj+ncq)*npv,nfe*npf*ncht,neNLb);
    
    % Compute global quantities to assemble element by element
    for ie=1:neNLb
        
        AA_K = AA(:,:,ie);
        
        zv = zeros(npv);
        zvf = zeros(npv,npf*nfe);
        zv3 = [zv zv zv];
        zv3f = zeros(3*npv,npf*nfe);
        zf2 = zeros(npf*nfe*nch,npf*nfe);
        
        
        Y1 = [AAmu(:,:,ie,1) zv zv zv zv zv; zv AAmu(:,:,ie,2) zv zv zv zv; zv zv AAmu(:,:,ie,3) zv zv zv;...
            zv zv zv alpha*AA_K zv zv; zv zv zv zv alpha*AA_K zv; zv zv zv zv zv alpha*AA_K];
        
        F_K = [AA_K zv zv; zv AA_K zv; zv zv AA_K];
        
        B_K = [zv BB(:,:,ie,3) -BB(:,:,ie,2); ...
            -BB(:,:,ie,3) zv BB(:,:,ie,1); ...
            BB(:,:,ie,2) -BB(:,:,ie,1) zv];
        
        P_K = [BB(:,:,ie,1); BB(:,:,ie,2) ; BB(:,:,ie,3)];
        
        Y2 = [-B_K zv3'; -delta*F_K -beta2*P_K];
        Y3 = [B_K.' -theta*F_K; zv3 -P_K.'];
        
        
        
        D_K = [DD(:,:,ie,1) -DD(:,:,ie,4) -DD(:,:,ie,5); ...
            -DD(:,:,ie,4)  DD(:,:,ie,2) -DD(:,:,ie,6); ...
            -DD(:,:,ie,5) -DD(:,:,ie,6)  DD(:,:,ie,3)];
        
        D_K = tau_t*D_K - omega^2*[AAep(:,:,ie,1) zv zv; ...
            zv AAep(:,:,ie,2) zv; ...
            zv zv AAep(:,:,ie,3)];
        H_K = AA_K + tau_n*HH(:,:,ie);
        Y4 = [D_K zv3'; zv3 H_K];
        
        
        C_K = [CC(:,:,ie,1) CC(:,:,ie,4); ...
            CC(:,:,ie,2) CC(:,:,ie,5); ...
            CC(:,:,ie,3) CC(:,:,ie,6)];
        
        G_K = [GG(:,:,ie,1); GG(:,:,ie,2); GG(:,:,ie,3)];
        Y5 = [-C_K zv3f;zv3f zv3f beta2*G_K];
        
        E_K = tau_t*[EE(:,:,ie,1) EE(:,:,ie,4); ...
            EE(:,:,ie,2) EE(:,:,ie,5); ...
            EE(:,:,ie,3) EE(:,:,ie,6)];
        
        N_K = tau_n*NN(:,:,ie);
        
        Y6 = [-E_K zv3f;zvf zvf -N_K];
        
        L_K = tau_t*[LL(:,:,ie,1) LL(:,:,ie,2) LL(:,:,ie,3); ...
            LL(:,:,ie,4) LL(:,:,ie,5) LL(:,:,ie,6)];
        
        R_K = [RR(:,:,ie,1) RR(:,:,ie,2) RR(:,:,ie,3); ...
            RR(:,:,ie,4) RR(:,:,ie,5) RR(:,:,ie,6)];
        
        Y7 = [-R_K zeros(npf*nfe*nch,3*npv); zeros(npf*nfe,3*npv) G_K.'];
        Y8 = [-L_K zeros(npf*nfe*nch,npv); zeros(npf*nfe,3*npv) -N_K.'];
        
        T_K = tau_n*TT(:,:,ie);
        M_K = tau_t*[MM(:,:,ie,1) MM(:,:,ie,2); ...
            MM(:,:,ie,2) MM(:,:,ie,3)];
        
        Y9 = [M_K zf2; zf2' T_K];
        
        
        YY = [Y1 Y2;Y3 Y4];
        tm = -YY\[Y5; Y6];
        zNL(:,:,ie) = tm;
        wNL(:,:,ie) = Y9 + [Y7 Y8]*tm;
    end
    ZNL(:,:,e1:e2) = zNL;
    WNL(:,:,e1:e2) = wNL;
    
    clear AA AAmu AAep BB CC DD EE LL RR MM TT GG HH zNL wNL
end

%% Right-hand-side
fprintf('rhs...')

SL = zeros(nfe*npf*nch,neL);
SNL = zeros(nfe*npf*ncht,neNL);
dirU = zeros(nfe*npf*ncht,neNL);

boundaryElementFace = cell2mat(mesh.boundaryElementFace);
t = 1;
for ie = boundaryElementFace(1,:)
    jf = boundaryElementFace(2,t);
    J = ((jf-1)*npf+1):jf*npf;
    I = perm(:,jf);
    
    nx = Nx(:,ie,jf);
    ny = Ny(:,ie,jf);
    nz = Nz(:,ie,jf);
    ef = shapft*ep(I,ie,1);
    tmp_ = shapfc*diag(gwfc.*detJf(:,ie,jf));
    
    ib = bcm_t(-bf(jf,ie));
    if ib > 2
        omegaeff = omega*sqrt(ef);
    else
        omegaeff = omega;
    end
    
    [islocal,qq] = ismember(ie,eleL);
    
    if islocal
        
        [fb1,fb2] = fbou_3d(ib,[nx ny nz],pfg(:,:,ie,jf),omegaeff,setup);
        SL(J,qq) = SL(J,qq) + tmp_*fb1;
        SL(J+npf*nfe,qq) = SL(J+npf*nfe,qq) + tmp_*fb2;
        
    else
        [~,qq] = ismember(ie,eleNL);
        
        [fb1,fb2] = fbou_3d(ib,[nx ny nz],pfg(:,:,ie,jf),omegaeff,setup);
        SNL(J,qq) = SNL(J,qq) + tmp_*fb1;
        SNL(J+npf*nfe,qq) = SNL(J+npf*nfe,qq) + tmp_*fb2;
        
        ib_n = bcm_n(-bf(jf,ie));
        if ib_n == 1           
            dirU(J+2*npf*nfe,qq) = 1;
        end
        
    end
    
    t = t + 1;
end

WL = reshape(WL,[nfe*npf,nch,nfe*npf,nch,neL]);
WL = permute(WL,[2 1 4 3 5]);
WNL = reshape(WNL,[nfe*npf,ncht,nfe*npf,ncht,neNL]);
WNL = permute(WNL,[2 1 4 3 5]);
SL = reshape(SL,[nfe*npf nch neL]);
SL = permute(SL,[2 1 3]);
SNL = reshape(SNL,[nfe*npf ncht neNL]);
SNL = permute(SNL,[2 1 3]);

dirU = reshape(dirU,[nfe*npf ncht neNL]);
dirU = permute(dirU,[2 1 3]);

%% Assembly of linear system using only dofs on the faces
fprintf('and assembly...')

elcon = reshape(mesh.elcon,npf*nfe,ne);
if neNL == 0
    
    elconL = zeros(npf*nfe*nch,ne);
    iL = zeros(nfe*npf*nch,nfe*npf*nch,neL);
    jL = zeros(nfe*npf*nch,nfe*npf*nch,neL);
    t = 1;
    for ie=eleL
        con = repmat((elcon(:,ie)'-1)*nch,nch,1)+repmat((1:nch)',1,nfe*npf);
        con = reshape(con,nfe*npf*nch,1);
        elconL(:,ie) = con;
        iL(:,:,t) = repmat(con ,1,nfe*npf*nch);
        jL(:,:,t) = repmat(con',nfe*npf*nch,1);
        t = t+1;
    end
    
    ndof = nch*npf*nf;
    WW = sparse(iL(:),jL(:),WL(:),ndof,ndof);
    tmL = iL(:,1,:);
    SS = sparse(tmL(:),1,SL(:),ndof,1);
    dof_solve = logical(sparse(tmL(:),1,ones(nfe*npf*nch*neL,1),ndof,1));
    
else
    

    interFace = mesh.interFace;
    faceNL = mesh.faceNL;
    faceL = mesh.faceL;
    
    ndofNL = npf*(numel(faceNL) + numel(interFace));
    ndofL = npf*numel(faceL);
    ndof = ncht*ndofNL + nch*ndofL;
    
    % just elements at the interface, treat carefully
    eleBou = mesh.f(interFace,end-1:end);
    nbouele = size(eleBou,1);
    elconL = zeros(npf*nfe*nch,ne); 
    eleBouL = zeros(1,nbouele);
    for i = 1:nbouele
        
        e1 = eleBou(i,1);
        e2 = eleBou(i,2);
        if mesh.metal(e2)
            e1 = eleBou(i,2);
            e2 = eleBou(i,1);
        end
        eleBouL(i) = e2;
        [~,~,tm2]=intersect(elcon(:,e1),elcon(:,e2));
        
        con = repmat((elcon(:,e2)'-1)*nch,nch,1)+repmat((1:nch)',1,nfe*npf);
        con = reshape(con,nfe*npf*nch,1);
        elconL(:,e2) = con;
        
        con_ = repmat((tm2'-1)*nch,nch,1)+repmat((1:nch)',1,npf);
        con_ = reshape(con_,npf*nch,1);
        
        new_ind = bsxfun(@minus,ncht*elcon(tm2,e2),nch:-1:1)';
        elconL(con_,e2) = new_ind(:);
        
        tm1_ = setdiff(1:nch*nfe*npf,con_);
        elconL(tm1_,e2) = elconL(tm1_,e2) + ndofNL*(ncht-nch);
    end
    for ie = setdiff(eleL,eleBouL)
        con = repmat((elcon(:,ie)'-1)*nch,nch,1)+repmat((1:nch)',1,nfe*npf);
        con = reshape(con,nfe*npf*nch,1);
        elconL(:,ie) = con + ndofNL*(ncht-nch);
    end
    
    % assembly indices for local and nonlocal faces
    iL = zeros(nfe*npf*nch,nfe*npf*nch,neL);
    jL = zeros(nfe*npf*nch,nfe*npf*nch,neL);
    t = 1;
    for ie=eleL
        con = elconL(:,ie);
        iL(:,:,t) = repmat(con ,1,nfe*npf*nch);
        jL(:,:,t) = repmat(con',nfe*npf*nch,1);
        t = t+1;
    end
    
    iNL = zeros(nfe*npf*ncht,nfe*npf*ncht,neNL);
    jNL = zeros(nfe*npf*ncht,nfe*npf*ncht,neNL);
    t = 1;
    elconNL = zeros(npf*nfe*ncht,ne);
    for ie=eleNL
        con = repmat((elcon(:,ie)'-1)*ncht,ncht,1)+repmat((1:ncht)',1,nfe*npf);
        con = reshape(con,nfe*npf*ncht,1);
        elconNL(:,ie) = con;
        iNL(:,:,t) = repmat(con ,1,nfe*npf*ncht);
        jNL(:,:,t) = repmat(con',nfe*npf*ncht,1);
        t = t +1;
    end
    
    % assemble large sparse system
    WW = sparse(iL(:),jL(:),WL(:),ndof,ndof) + ...
        sparse(iNL(:),jNL(:),WNL(:),ndof,ndof);
    tmL = iL(:,1,:);
    tmNL = iNL(:,1,:);
    SS = sparse(tmL(:),1,SL(:),ndof,1) + sparse(tmNL(:),1,SNL(:),ndof,1);
    
    dirUU = logical(sparse(tmNL(:),1,dirU(:),ndof,1));
    SS(dirUU) = [];
    WW(dirUU,:) = [];
    WW(:,dirUU) = [];
    dof_solve = logical(1-dirUU);
end
fprintf('done in %.1f seconds\nSolving linear system...',toc)
tic

clear WL WNL SL SNL

%% Solution of linear system

UH = zeros(ndof,1);
UH(dof_solve) = WW\SS;

fprintf('done in %.1f seconds\n',toc)

clear WW SS


%% Reconstruction of local variables from global variables
fprintf('Recovering local variables...')

EDG = zeros(npv,ne,ncu);
VDG = zeros(npv,ne,ncv);
JDG = zeros(npv,ne,ncj);
QDG = zeros(npv,ne);

for i = 1:neNL
    ie = eleNL(i);
    imap = elconNL(:,ie);
    
    uh = permute(reshape(UH(imap),ncht,[]),[2 1]);
    
    loc = ZNL(:,:,i)*uh(:);
    
    vdg = loc(1:ncv*npv);
    jdg = loc(ncv*npv+1:(ncv+ncj)*npv);
    udg = loc((ncv+ncj)*npv+1:(ncv+ncj+ncu)*npv);
    qdg = loc((ncv+ncj+ncu)*npv+1:end);
    
    
    VDG(:,ie,1) = vdg(1:npv);
    VDG(:,ie,2) = vdg(npv+1:2*npv);
    VDG(:,ie,3) = vdg(2*npv+1:3*npv);
    EDG(:,ie,1) = udg(1:npv);
    EDG(:,ie,2) = udg(npv+1:2*npv);
    EDG(:,ie,3) = udg(2*npv+1:3*npv);
    QDG(:,ie) = qdg;
    JDG(:,ie,1) = jdg(1:npv,:);
    JDG(:,ie,2) = jdg(npv+1:2*npv,:);
    JDG(:,ie,3) = jdg(2*npv+1:3*npv,:);
    
    
end


for i = 1:neL
    ie = eleL(i);
    imap = elconL(:,ie);
    
    uh = permute(reshape(UH(imap),nch,[]),[2 1]);
    
    loc = ZL(:,:,i)*uh(:);
    
    vdg = loc(1:ncv*npv);
    udg = loc(ncv*npv+1:(ncv+ncu)*npv);
    
    
    VDG(:,ie,1) = vdg(1:npv);
    VDG(:,ie,2) = vdg(npv+1:2*npv);
    VDG(:,ie,3) = vdg(2*npv+1:3*npv);
    EDG(:,ie,1) = udg(1:npv);
    EDG(:,ie,2) = udg(npv+1:2*npv);
    EDG(:,ie,3) = udg(2*npv+1:3*npv);
    
    
    
end

EDG = permute(EDG,[1 3 2]);
HDG = -1i/omega*permute(VDG,[1 3 2]);
JDG = permute(JDG,[1 3 2]);

if ishydro    
    RDG = -1i/omega*permute(QDG,[1 3 2]);
else
    in = mesh.metal == 1;
    JDG(:,:,in) = 1i*omegap^2/(omega+1i*gamma)*EDG(:,:,in);
    RDG = [];
end

fprintf('done!')

end


function [fb1,fb2] = fbou_3d(ib,nl,pg,omega,data)
%FHAT flux function

%ng = size(pg,1);
x  = pg(:,1);
y  = pg(:,2);
z  = pg(:,3);
nx = nl(:,1);
ny = nl(:,2);
nz = nl(:,3);
[t1x,t1y,t1z,t2x,t2y,t2z] = tangentvectors(nx,ny,nz);
k = data.k;
p = data.p;

switch ib
     case 1  % PEC                
        fb1 = 0*x;                
        fb2 = 0*x;           
    case 2  % PMC                 
        fb1 = 0*x;                
        fb2 = 0*x;   

    case {3,4}  % Silver-Muller:  curl E cross n - i omega E^t  =  curl E^in cross n - i omega n cross E^in cross n
        

        kxp = cross(k,p);
        
        tm1 = kxp(1).*(t1z.*ny-t1y.*nz) + kxp(2).*(t1x.*nz-t1z.*nx) + kxp(3).*(t1y.*nx-t1x.*ny) - (p(1).*t1x + p(2).*t1y + p(3).*t1z);
        tm2 = kxp(1).*(t2z.*ny-t2y.*nz) + kxp(2).*(t2x.*nz-t2z.*nx) + kxp(3).*(t2y.*nx-t2x.*ny) - (p(1).*t2x + p(2).*t2y + p(3).*t2z);
        
        
        e  = 1i*omega.*exp(1i*omega.*(k(1)*x+k(2)*y+k(3)*z));
        
        fb1 = e.*tm1;
        fb2 = e.*tm2;
        
    
    otherwise
        error('unknown boundary type');
end
end

function [t1x,t1y,t1z,t2x,t2y,t2z] = tangentvectors(nx,ny,nz)
% determine two tangent vectors on a face from the normal vector


ngf = size(nx,1);
nm = [mean(abs(nx)); mean(abs(ny)); mean(abs(nz))];
[~,id] = max(nm+repmat([1e-12;2e-12;3e-12],[1 size(nx,2)])); 
t1x = nx;
t1y = ny;
t1z = nz;
t2x = nx;
t2y = ny;
t2z = nz;

i1 = find(id==1);
t1x(:,i1) = -ny(:,i1)./nx(:,i1);
t1y(:,i1) = ones(ngf,length(i1));
t1z(:,i1) = zeros(ngf,length(i1));
t2x(:,i1) = -nz(:,i1)./nx(:,i1);
t2y(:,i1) = zeros(ngf,length(i1));
t2z(:,i1) = ones(ngf,length(i1));

i1 = find(id==2);
t1x(:,i1) = ones(ngf,length(i1));
t1y(:,i1) = -nx(:,i1)./ny(:,i1);
t1z(:,i1) = zeros(ngf,length(i1));
t2x(:,i1) = zeros(ngf,length(i1));
t2y(:,i1) = -nz(:,i1)./ny(:,i1);
t2z(:,i1) = ones(ngf,length(i1));

i1 = find(id==3);
t1x(:,i1) = ones(ngf,length(i1));
t1y(:,i1) = zeros(ngf,length(i1));
t1z(:,i1) = -nx(:,i1)./nz(:,i1);
t2x(:,i1) = zeros(ngf,length(i1));
t2y(:,i1) = ones(ngf,length(i1));
t2z(:,i1) = -ny(:,i1)./nz(:,i1);    

end

