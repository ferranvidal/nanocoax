function h = faceplot(mesh,UDG,zval,plotmesh,clim)

% plotmesh is whether the 2D is plotted (1) or not(0)
perm = mesh.perm;
ne = mesh.ne;
nf = mesh.nf;
nd = mesh.nd;
[npf,nfe] = size(perm);
elcon = reshape(mesh.elcon,[npf nfe ne]);

in = mesh.f(:,end)<0;
ns = sum(in);
dgnodes = zeros(npf,nd,ns);
u = zeros(npf,ns);
if isempty(UDG)
    UDG = zeros(mesh.npv,ne);
    ismesh = 1;
else
    UDG = reshape(UDG,[],ne);
    ismesh = 0;
end

z = unique(mesh.p(:,3));
[~,i]= min(abs(z-zval));
fprintf('Showing solution at z = %.4f\n',z(i))
func = @(p) all(abs(p(:,3)-z(i))<1e-5);


j = 0;
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i
    tf = mesh.f(i,1:end-2);
    pf = (mesh.p(tf,:));
    in = feval(func,pf);
    if in == 1
        %         hold on
        j = j + 1;
        if fi(2)>0
            kf = mesh.t2f(fi(1),:);  % obtain neighboring faces
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;
            kf = mesh.t2f(fi(2),:);    % obtain neighboring faces
            i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element
            j2 = elcon(:,i2,fi(2)) - (i-1)*npf;
            u1 = UDG(perm(j1,i1),fi(1));
            u2 = UDG(perm(j2,i2),fi(2));
            u(:,j) = (u1+u2)/2;
            dgnodes(:,:,j) = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));
            %[mesh.dgnodes(perm(j1,i1),1:nd,fi(1))-mesh.dgnodes(perm(j2,i2),1:nd,fi(2))]
            %pause
        else
            kf = mesh.t2f(fi(1),:);    % obtain neighboring faces
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;
            u(:,j) = UDG(perm(j1,i1),fi(1));
            dgnodes(:,:,j) = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));
        end
        
        
        %         plot3(pf(:,1),pf(:,2),pf(:,3),'r*')
        %         plot3(dgnodes(:,1,j),dgnodes(:,2,j),dgnodes(:,3,j),'b.')
    end
end

dgnodes = dgnodes(:,:,1:j);
u = u(:,1:j);

cmapbool = all(u(:) > 0);

plocal=mesh.plocfc;
tlocal=mesh.tlocfc;

dgnodes = permute(dgnodes,[1,3,2]);
nt = size(u,2);
npln=size(plocal,1);
nodesvis=reshape(dgnodes,[npln*nt,nd]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

nodesvis(:,3) = [];

if ~ ismesh
    h = patch('vertices',nodesvis,'faces',tvis,'cdata',u(:), ...
        'facecol','interp','edgec','none');
end
if cmapbool
    colormap inferno
    edgec = [.8 .8 .8];
else
    colormap bluewhitered
    edgec = 'k';
end

if (nargin >= 4 && plotmesh==1) || ismesh
    pars={'facecolor','none','edgecolor',edgec,'Linew',0.5};
    e=boundedges(plocal,tlocal,mesh.elemtype);
    e1=segcollect(e);
    axis equal,axis off
    hh=zeros(nt,1);
    for it=1:nt
        px=dgnodes(:,it,1);
        py=dgnodes(:,it,2);
        
        pz=0*px;
        hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0.0*e1{1}',pars{:});
    end
end
if ~ ismesh
    colorbar
    axis off
    if nargin>=5 && ~isempty(clim)
        set(gca,'clim',clim);
    end
end
set(gca,'fontsize',16)
set(gcf,'color','w')


function e=boundedges(p,t,elemtype)
%BOUNDEDGES Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T)

% Form all edges, non-duplicates are boundary edges

if elemtype==0
    edges=[t(:,[1,2]);
        t(:,[1,3]);
        t(:,[2,3])];
    node3=[t(:,3);t(:,2);t(:,1)];
else
    edges=[t(:,[1,2]);
        t(:,[2,3]);
        t(:,[3,4]);
        t(:,[4,1]);];
    node3=[t(:,4);t(:,3);t(:,2);t(:,1)];
end
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);
node3=node3(ix(qx));

% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(ix,[1,2])=e(ix,[2,1]);


function e1=segcollect(e)
%SEGCOLLECT Collect polygons from edge segments.

ue=unique(e(:));
he=histc(e(:),ue);
current=ue(min(find(he==1))); % Find an endpoint
if isempty(current) % Closed curve
    current=e(1,1);
end
e1=current;
while ~isempty(e)
    ix=min(find(e(:,1)==e1(end)));
    if isempty(ix)
        ix=min(find(e(:,2)==e1(end)));
        if isempty(ix) % >1 disjoint curves, recur
            rest=segcollect(e);
            e1={e1,rest{:}};
            return;
        end
        next=e(ix,1);
    else
        next=e(ix,2);
    end
    e1=[e1,next];
    e(ix,:)=[];
end
e1={e1};

