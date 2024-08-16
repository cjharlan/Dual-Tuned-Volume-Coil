function [B]=BField2(R1p,R2p,Ro,tol);

echo on
%This function computes the Biot-Savart integral over a 
% length of wire extending from R1p to R2p, to find the
% magnetic field intensity at Ro, the observation point.
%
%Usage: [B]=BField(R1p,R2p,Ro);
%
%
% R1p    = [x, y, z] of POSITIVE edge of wire (current source)
% R2p    = [x, y, z] of NEGATIVE edge of wire (current sink)
% Ro     = [x, y, z] coordinates of the observation point.
%          [n, NPE, NOP] n=1..3 for x,y,z positions of imaging plane
% B      = [Bx, By, Bz] at the observation point(s). (A/m^2)
%
% ALL LOCATION DIMENSIONS IN METERS
% 
echo off

if (nargin<4)
    tol=1e-3;
end

mu=4*pi*1e-7;
[sxo syo szo] = size(Ro);
dl = [R2p(1)-R1p(1) R2p(2)-R1p(2) R2p(3)-R1p(3)];
clength = norm(dl);
dl = dl/clength;
intfac=256;

if sxo==szo
    rxp = linspace(R1p(1),R2p(1),intfac);
    ryp = linspace(R1p(2),R2p(2),intfac);
    rzp = linspace(R1p(3),R2p(3),intfac);
    rx  = Ro(1)-rxp;
    ry  = Ro(2)-ryp;
    rz  = Ro(3)-rzp;
    R   = sqrt(rx.^2+ry.^2+rz.^2);
    Bx  = sum(((dl(2)*rz-dl(3)*ry)./(R.^3)).*(R>tol));
    By  = sum(((dl(3)*rx-dl(1)*rz)./(R.^3)).*(R>tol));
    Bz  = sum(((dl(1)*ry-dl(2)*rx)./(R.^3)).*(R>tol));
    B   = [Bx By Bz]*mu/4/pi*clength/intfac;
else
    rxp = linspace(R1p(1),R2p(1),intfac);
    ryp = linspace(R1p(2),R2p(2),intfac);
    rzp = linspace(R1p(3),R2p(3),intfac);
    Bx=zeros(szo,syo);
    By=Bx;
    Bz=Bx;
    for ii=1:intfac
        rx = squeeze(Ro(1,:,:))-rxp(ii);
        ry = squeeze(Ro(2,:,:))-ryp(ii);
        rz = squeeze(Ro(3,:,:))-rzp(ii);
        R   = sqrt(rx.^2+ry.^2+rz.^2);
        Bx  = Bx + ((dl(2)*rz-dl(3)*ry)./(R.^3)).*(R>tol);
        By  = By + ((dl(3)*rx-dl(1)*rz)./(R.^3)).*(R>tol);
        Bz  = Bz + ((dl(1)*ry-dl(2)*rx)./(R.^3)).*(R>tol);
    end
    B=zeros(3,syo,szo);
    B(1,:,:)=Bx;
    B(2,:,:)=By;
    B(3,:,:)=Bz;
    B=B*mu/4/pi*clength/intfac;
end
