function [ax, sag, cor, axl, axr] = birdfield(blen, bdiam, nleg, npts, fovz, fovx);

l1=bdiam*0.6;
l2=blen*0.6;
tol=5e-4;

if nargin<4
    npts=65;
else
    if mod(npts,2)==0
        npts=npts+1;
    end
end

axl=linspace(l2,-l2,npts);
axr=linspace(-l1,l1,npts);
afov=zeros(3,npts,npts);
afov(2,:,:)=(linspace(l1,-l1,npts).')*ones(1,npts); % y-dir
afov(1,:,:)=ones(npts,1)*linspace(-l1,l1,npts);         % x-dir
sfov=zeros(3,npts,npts);
sfov(3,:,:)=(linspace(l2,-l2,npts).')*ones(1,npts); % z-dir
sfov(2,:,:)=ones(npts,1)*linspace(-l1,l1,npts);         % y-dir
cfov=zeros(3,npts,npts);
cfov(3,:,:)=(linspace(l2,-l2,npts).')*ones(1,npts); % z-dir
cfov(1,:,:)=ones(npts,1)*linspace(-l1,l1,npts);         % x-dir

ax=zeros(npts);
sag=zeros(npts);
cor=zeros(npts);

zind=ceil(npts/2);
[tmp zmin]=min(abs(axl-fovz/2));
[tmp zmax]=min(abs(axl+fovz/2));
[tmp xmin]=min(abs(axr+fovz/2));
[tmp xmax]=min(abs(axr-fovz/2));
zmask=zeros(1,npts);
zmask(zmin:zmax)=ones(1,zmax-zmin+1);
xmask=zeros(1,npts);
xmask(xmin:xmax)=ones(1,xmax-xmin+1);

fprintf('Segment:');
for ii=0:nleg-1
    %ileg = cos(2.0*pi*ii/nleg)-cos(2.0*pi*(ii-1)/nleg);
    %iend = cos(2.0*pi*ii/nleg);
    %Shit!
    %ileg = cos(2.0*pi*ii/nleg);
    ileg = cos(2.0*pi*ii/nleg)-cos(2.0*pi*(ii-1)/nleg)+...
        (-j)*sin(2.0*pi*ii/nleg)-(-j)*sin(2.0*pi*(ii-1)/nleg);
    iend = cos(2.0*pi*ii/nleg) + (-j)*sin(2.0*pi*ii/nleg);
    %iend=1;
    %ileg=0;
    %fprintf('%02d.. Ileg = %5.3f; Iend = %5.3f\n',ii+1,ileg,iend);
    fprintf('%02d.. ',ii+1);
    if abs(ileg)>1e-10
        legx = .5*bdiam*cos(2*pi*ii/nleg);
        legy = .5*bdiam*sin(2*pi*ii/nleg);
        [Bt] = bfield2([legx, legy, -blen/2],[legx, legy, blen/2],afov,tol);
        ax=ax+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*ileg;
        [Bt] = bfield2([legx, legy, -blen/2],[legx, legy, blen/2],sfov,tol);
        sag=sag+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*ileg;
        [Bt] = bfield2([legx, legy, -blen/2],[legx, legy, blen/2],cfov,tol);
        cor=cor+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*ileg;
    end
    if abs(iend)>1e-10
        legxs = .5*bdiam*cos(2*pi*ii/nleg);
        legys = .5*bdiam*sin(2*pi*ii/nleg);
        legxe = .5*bdiam*cos(2*pi*(ii+1)/nleg);
        legye = .5*bdiam*sin(2*pi*(ii+1)/nleg);
        [Bt]=bfield2([legxs,legys,blen/2],[legxe,legye,blen/2],afov,tol);
        ax=ax+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
        [Bt]=bfield2([legxs,legys,-blen/2],[legxe,legye,-blen/2],afov,tol);
        ax=ax-squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
        [Bt]=bfield2([legxs,legys,blen/2],[legxe,legye,blen/2],sfov,tol);
        sag=sag+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
        [Bt]=bfield2([legxs,legys,-blen/2],[legxe,legye,-blen/2],sfov,tol);
        sag=sag-squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
        [Bt]=bfield2([legxs,legys,blen/2],[legxe,legye,blen/2],cfov,tol);
        cor=cor+squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
        [Bt]=bfield2([legxs,legys,-blen/2],[legxe,legye,-blen/2],cfov,tol);
        cor=cor-squeeze(Bt(1,:,:)+j*Bt(2,:,:))*iend;
    end
end

figure(1)
imagesc([-l1,l1],[-l1,l1],abs(ax),[0 max(max(abs(ax)))]);
colorbar('vert')
title('Axial field plot')

figure(2)
contourf(axr,axr,abs(ax),10);
colorbar('vert')
title('Axial field plot')

figure(3)
imagesc([-l1,l1],[-l2,l2],abs(sag),[0 max(max(abs(sag)))]);
colorbar('vert')
title('Saggital field plot')
figure(4)
imagesc([-l1,l1],[-l2,l2],abs(cor),[0 max(max(abs(cor)))]);
colorbar('vert')
title('Coronal field plot')

figure(5)
subplot(211)
plot(axr,squeeze(abs(ax(zind,:))));
axis([axr(1) axr(end) 0 max(squeeze(abs(ax(zind,:))))]);
title('Sensitivity along x');
subplot(212)
plot(axl,squeeze(abs(sag(:,zind))));
axis([axl(end) axl(1) 0 max(squeeze(abs(sag(:,zind))))]);
title('Sensitivity along z');

figure(6)
subplot(211)
refval=ax(zind,zind);
unifxp=abs(ax(:,zind)-refval)/abs(refval);
plot(100*axr,unifxp,100*axr,unifxp.*(xmask.'));
title('Percent deviation from center point');
ylabel('Along X-Axis');
axis([axr(1)*100 axr(end)*100 0 0.2]);
subplot(212)
unifzp=abs(sag(:,zind)-refval)/abs(refval);
plot(100*axl,unifzp,100*axl,unifzp.*(zmask.'));
ylabel('Along Z-Axis');
xlabel('(cm)');
grid on
axis([axl(end)*100 axl(1)*100 0 0.2]);

unifz = max(abs(sag(zmin:zmax,zind)-refval))/abs(refval);
unifx = max(abs(ax(xmin:xmax,zind)-refval))/abs(refval);
unify = max(abs(ax(zind,xmin:xmax)-refval))/abs(refval);
fprintf('\nXDev: %5.3f (%3.1fcm); YDev: %5.3f (%3.1fcm); ZDev: %5.3f (%3.1fcm\n)',...
    unifx,abs(axr(xmax)-axr(xmin))*100,...
    unify,abs(axr(xmax)-axr(xmin))*100,...
    unifz,abs(axl(zmax)-axl(zmin))*100);
