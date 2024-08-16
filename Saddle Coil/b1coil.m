
function [Bt]=b1coil(xo,yo,theta,basisr);

echo on
%This function computes the Biot-Savart integral for a circular coil,
%radius r, offset from the center of a field of view by xo, yo, with normal
%at angle theta from the vertical axis of the FOV.
% on the surface of a cylinder.  
%
%Usage: [B]=b1coil(xo,yo,theta,basisr);
%
%
% B      = [Bx, By, Bz] at the observation point(s). (A/m^2)
% theta  = theta angle from the vertical axis of the FOV
% basisr = [n, m] m=1..3 for x,y,z coords of n obs points
%
% ** ALL DIMENSIONS ARE METERS AND RADIANS! **
% 
echo off

mu=4*pi*1e-7;
tol=1e-4; % do not calculate less than 0.1mm from conductor.
coilr=.010/2; %10-mm diameter
intfac=32;
clength=2*pi*coilr/intfac;
nbf=length(basisr);% How many basis functions to calc?
Bt=zeros(1,nbf);



%Calculate distance from center of FOV to each obs point, given angle:
Ro=zeros(size(basisr));
basisr(:,1)=basisr(:,1)-xo;
basisr(:,2)=basisr(:,2)-yo;
%Rotate:
Ro=basisr(:,1)*cos(-theta)-basisr(:,2)*sin(-theta);
Ro(:,2)=basisr(:,1)*sin(-theta)+basisr(:,2)*cos(-theta);
Ro(:,3)=0;
%Leave rotation along z alone...for now.
%%Now add offset to coil:
%Ro(:,1)=Ro(:,1)-xo;
%Ro(:,2)=Ro(:,2)-yo;
%Ro(:,3)=0; %No offset along Z...for now

for ii=1:intfac
    % Calculate position and orientation of conductors, relative to
    % the center of the coil
    thp = 2*pi*(ii-1)/intfac;
    Rp=coilr*([cos(thp) 0 sin(thp)]);
    dlx=-sin(thp);
    dly=0;
    dlz=cos(thp);
    % Calculate the contribution of this segment to obs points:
    for jj=1:nbf
        rx  = Ro(jj,1)-Rp(1);
        ry  = Ro(jj,2)-Rp(2);
        rz  = Ro(jj,3)-Rp(3);
        R   = sqrt(rx.^2+ry.^2+rz.^2);
        Bx  = sum(((dly.*rz-dlz.*ry)./(R.^3)).*(abs(R)>tol));
        By  = sum(((dlz.*rx-dlx.*rz)./(R.^3)).*(abs(R)>tol));
        %Bz  = sum(((dlx.*ry-dly.*rx)./(R.^3)).*(abs(R)>tol));
        Bt(jj)   = Bt(jj) + (Bx +1i*By);
    end
end
Bt=Bt*mu/4/pi*clength/intfac;