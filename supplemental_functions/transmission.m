function [S,A] = transmission(mesh,EDG,HDG,func)

% Computes \int_A \abs{Re(ExH*).n}
if isempty(HDG), S = 0; A = 0; return; end

perm = mesh.perm;
ne = mesh.ne;
nf = mesh.nf;
nd = mesh.nd;
[npf,nfe] = size(perm);
elcon = reshape(mesh.elcon,[npf nfe ne]);


if size(EDG,2) > nd
    EDG = permute(EDG,[1 3 2]) ;
end
if size(HDG,2) > nd
    HDG = permute(HDG,[1 3 2]) ;
end


dgnodes = mesh.dgnodes;


in = mesh.f(:,end)<0;
ns = sum(in);
X = zeros(npf,nd,ns);

E = zeros(npf,nd,ns);
H = zeros(npf,nd,ns);

j = 0;
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i
    tf = mesh.f(i,1:end-2);
    pf = (mesh.p(tf,:));
    in = feval(func,pf);
    
    
    if (in == 1)
        j = j + 1;
        if fi(2)>0
            kf = mesh.t2f(fi(1),:);  % obtain neighboring faces
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;
            kf = mesh.t2f(fi(2),:);    % obtain neighboring faces
            i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element
            j2 = elcon(:,i2,fi(2)) - (i-1)*npf;
            u1 = EDG(perm(j1,i1),:,fi(1));
            u2 = EDG(perm(j2,i2),:,fi(2));
            E(:,:,j) = (u1+u2)/2;
            u1 = HDG(perm(j1,i1),:,fi(1));
            u2 = HDG(perm(j2,i2),:,fi(2));
            H(:,:,j) = (u1+u2)/2;
            X(:,:,j) = dgnodes(perm(j1,i1),1:nd,fi(1));

        else
            kf = mesh.t2f(fi(1),:);    % obtain neighboring faces
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;
            E(:,:,j) = EDG(perm(j1,i1),:,fi(1));
            H(:,:,j) = HDG(perm(j1,i1),:,fi(1));
            X(:,:,j) = dgnodes(perm(j1,i1),1:nd,fi(1));
        end
        
    end
end
X = X(:,:,1:j);
E = E(:,:,1:j);
H = H(:,:,1:j);
ns = j;

master = mesh.master;
gwfc = master.gwfc;
shapfc = master.shapfc;
[npf ,ngf, nd] = size(shapfc);
shapft  = shapfc(:,:,1)';
dshapft = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[ngf*(nd-1) npf]);


[ngf,~] = size(shapft);

A = 0;
S = 0;
for jf=1:ns
    pb = X(:,:,jf);
    dpg = permute(reshape(dshapft*pb,[ngf nd-1 nd]),[1 3 2]);
    nx = reshape(dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2),[ngf 1]);
    ny = reshape(dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2),[ngf 1]);
    nz = reshape(dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2),[ngf 1]);
    detf = sqrt(nx.^2+ny.^2+nz.^2);
    nx = nx./detf;
    ny = ny./detf;
    nz = nz./detf;
    nn = [nx ny nz];
    
    Eg = reshape(shapft*E(:,:,jf),[ngf nd]);
    Hg = reshape(shapft*H(:,:,jf),[ngf nd]);
    P = real(cross(Eg,conj(Hg),2));
    Pn = dot(P,nn,2);
    S  = S + sum(abs(Pn).*gwfc.*detf);
    A  = A + sum(gwfc.*detf);
end
