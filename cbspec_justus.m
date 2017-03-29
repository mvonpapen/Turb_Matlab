%% Numerically integrates P_i(f,\theta)=\int dk^3 E_i(k) \delta(2\pi f - k*v) for given kx(=2pi f/v, kx is in rotated frame with kx||v) and \theta under assumption of Taylor's hypothesis. Model of P_i(k) is critically balanced with k_para~L^(1/3)k_perp^(2/3) in inertial range and k_\para~L^(1/3)\rho^(1/3)k_\perp^(1/3) in kinetic range. This code is introduced and discussed in detail in: von Papen, 2014, PhD thesis, Cologne, Germany. More info on the topic can be found in Cho et al., 2002 (ApJ), Forman et al., 2011 (ApJ), or Wicks et al., 2012 (ApJ)
%%
%%
%%  Input:
%%      f       Frequncy in Hz
%%      theta   Field-to-flow angle in degrees
%%      n       Number of nodes in one dimension
%%      fun     function to evaluate ( 'exp', 'expdamp', 'gauss', 'gaussdamp', 'heavi' )
%%      L       Outer scale in m
%%      rho     Gyro radius (or else) in m
%%      v       Plasma bulk velocity in m/s
%%      va      Alfven velocity in m/s
%%      B       Magnetic field strength in nT
%%      bounds  Boundary in log10 for integration ( bounds = [kmin kmax] )
%%      ratioTP Ratio between toroidal T and poloidal P fluctuations
%%      rhoe    Electron gyro radius (or else)
%%      si_in   Spectral index in [MHD kinetic] range
%%      cb_in   Critical balance exponent in [MHD kinetic] range
%%
%%  Output:
%%      P       PSD in nT^2/Hz
%%
%%  Author: Michael von Papen
%%  Date: 04.08.2014
%%
%%  Notice, please acknowledge the use of this code in any publications.
%%  Reference:
%%      von Papen, M., 2014: Turbulence in Saturn's Magnetosphere and Forward Modeling of Reduced Spectra from Three-Dimensional Wave Vector Space, PhD thesis, Cologne
%%
%%  I would be pleased to receive a copy of such publications under:
%%      Michael von Papen
%%      Institute of Geophysics
%%      Pohligstr. 3
%%      50969 Cologne, Germany
%%      E-mail: vonpapen@geo.uni-koeln.de
%% -----------------------------------------------

function  [ P ] = cbspec ( f, theta, n, fun, L, rho, v, va, B, bounds, ratioTP, rhoe, si_in, cb_in)


%% Check Input
if nargin<14; cb_in=[2/3 1/3]; end
if nargin<13; si_in=[-10/3 -11/3]; end
if nargin<11; ratioTP = 0; end
if nargin<10; bounds=[-10 -2]; end
if nargin<9; B=1; end
if nargin<8; va=60e3; end
if nargin<7; v=600e3; end
if nargin<6; rho=1e5; end
if nargin<5; L=1e9; end
if nargin<4; fun='exp'; end
if nargin<3; n=500; end
if nargin<2; theta=[0 90]; end
if nargin<12; rhoe=rho/42.85; end


%% Wave vector in rotated coordinate system
kix=2*pi*f/v;

%% Basic Parameter
si  = si_in(1); %-10/3 -> k^{-5/3}
si2 = si_in(2); %-11/3 -> k^{-7/3}
cb  = cb_in(1); % 2/3  -> alfven
cb2 = cb_in(2); % 1/3  -> KAW


%% K-space gridpoints
kmin=bounds(1);
kmax=bounds(2);

% Range to cover in each dimension of the k-integration.
% n is the nunber of nodes integrated over for each sign in each direction,
% k'z covers both signs, there are thus 2n nodes in that direction.
kRange = 10.^(kmin + (0:n-1)*(kmax-kmin)/n);

% Log-distribute z-momenta over both positive and negative range.
kiz = repmat([-kRange(end:-1:1) kRange], n, 1);

% The y-component only needs positive values, due to symmetry.
% ky' is actually the same as ky, hence we don't use the name kiy here.
ky = repmat(transpose(kRange), 1, 2*n);



%% Ion-cyclotron frequency (wic) cut-off.
% All k_para ~ w_ic/V_A are subject to ion-cyclotron damping.
% Thus, parallel scales cannot reach k_para >> w_ic/V_A
wic = 1.6e-19*B*1e-9/1.67e-27; %wic for protons


%% Set output variable
if ratioTP ~=0
    P = zeros(length(kix),length(theta),4);
else
    P = zeros(length(kix),length(theta));
end

%% Begin with loop over theta
for k=1:length(theta)


  %% Begin loop over frequency
  for i=1:length(kix)

    %% Calculate PSD at z
    if theta(k)==0
      kx = kix(i)*sind(theta(k)) - kiz.*cosd(theta(k)); %=kx in unrotated system
      %ky=kiy in unrotated system
      kz = kix(i)*cosd(theta(k)) + kiz.*sind(theta(k)); %=kz in unrotated system
    else
                      %=kx in unrotated system with kiz=0 on x-axis
      kx = kix(i)*sind(theta(k))                            ...
          - kiz.*cosd(theta(k))                             ...
          + kix(i).*cosd(theta(k))^2/sind(theta(k));
      %ky=kiy in unrotated system
      kz = kiz.*sind(theta(k)); %=kz in unrotated system with kiz=0 on x-axis
    end

    kern = zeros(n, 2*n);

    kp2 = ky.^2 + kx.^2; % k_perp^2
    kabs2 = kp2 + kz.^2; % |k|

    %% Equations written in unprimed coordinates for the sake of brevity, but
    %  integration is done over primed variables, which is why dkiy and dkiz is used.

    %% Single components Alfven cascade.
    i1 = find(kp2 <= 1/rho^2 & kp2 > 1/L^2 & kabs2 <= 10^(2*kmax));
                                  % take out '& kp2 > 1/L^2' when checking for fmax
    
    switch fun
      case 'exp'
          kern(i1) = kp2(i1).^(si/2)                                                       ...
                         .* exp(-L^(1-cb) .* abs( kz(i1) ) ./ kp2(i1).^(cb/2) );
      case 'expdamp' % for Ti=Te => rho_e = sqrt(me/mi) * rho_i
          kern(i1) = kp2(i1).^(si/2)                                                       ...
                         .* exp( - L^(1-cb) .* abs(kz(i1) ./ kp2(i1).^(cb/2))              ...
                                 - sqrt( kp2(i1) ) * rhoe                                  ...
                                 - abs(kz(i1)) * va / wic );
      case 'gauss'
          kern(i1) = kp2(i1).^(si/2)                                                       ...
                         .* exp( - (L^(1-cb) * abs( kz(i1) ) ./ kp2(i1).^(cb/2) - 1 ).^2 ) ...
                         ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = kp2(i1).^(si/2)                                                       ...
                         .* exp( - ( L^(1-cb) * abs(kz(i1)) ./ kp2(i1).^(cb/2) - 1 ).^2    ...
                                 - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )         ...
                         ./ sqrt(pi); 
      case 'heavi'
          i2 = find( L^(1-cb) * abs(kz(i1) ./ kp2(i1).^(cb/2) ) <= 1 );
          kern( i1(i2) ) = kp2(i1(i2)).^(si/2);
    end


    %% Single components KAW cascade
    i1 = find(kp2 > 1/rho^2 & kp2 > 1/L^2 & kabs2 <= 10^(2*kmax) & kp2 <= 1/rhoe^2);
    switch fun
      case 'exp'
          kern(i1) = rho^(si2-si)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( -L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) );
      case 'expdamp' % for Ti=Te => rho_e=sqrt(me/mi)*rho_i
          kern(i1) = rho^(si2-si)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( - L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2)            ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic );
      case 'gauss'
          kern(i1) = rho^(si2-si)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( -(L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) - 1 ).^2 ) ...
                        ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = rho^(si2-si)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( - (L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) - 1 ).^2  ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )                      ...
                        ./ sqrt(pi);
      case 'heavi'
          i2 = find( L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) <= 1 );
          kern( i1(i2) ) = rho^(si2-si) .* kp2(i1(i2)).^(si2/2);
    end

    %% Single components of cascade at electron scales
    i1 = find(kp2 > 1/rhoe^2 & kabs2 <= 10^(2*kmax));
    switch fun
      case 'exp'
          kern(i1) = rho^(si2-si)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( -L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) * rhoe^cb2 );
      case 'expdamp' % for Ti=Te => rho_e=sqrt(me/mi)*rho_i
          kern(i1) = rho^(si2-si)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( - L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) * rhoe^cb2            ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic);
      case 'gauss'
          kern(i1) = rho^(si2-si)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( -(L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) * rhoe^cb2 - 1 ).^2 ) ...
                        ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = rho^(si2-si)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( - (L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) * rhoe^cb2 - 1 ).^2  ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )             ...
                        ./ sqrt(pi);
      case 'heavi'
          i2 = find( L^(1-cb) * rho^(cb-cb2) * abs(kz(i1)) * rhoe^cb2 <= 1 );
          kern( i1(i2) ) = rho^(si2-si) .* kp2(i1(i2)).^(si2/2);
    end

    
    % Algorithm for 2D-integration over the kxy-domain ranges.
    integral_dkxy = @(fk) trapz( ky(:,1), trapz(kiz(1,:), fk, 2) );

    %% Full version with Toroidal and Poloidal parts
    % Sum up to get power for one ky value
    if ratioTP ~= 0
      Tor = ratioTP / (1+ratioTP) * kern ./ kp2;
      
      Pol = kern / ratioTP ./ kabs2;

      %% Add Spectra to PSD
      P(i,k,1) = integral_dkxy( ky.^2.*Tor + (kx.*kz).^2 ./ kp2 .* Pol );
      P(i,k,2) = integral_dkxy( kx.^2.*Tor + (ky.*kz).^2 ./ kp2 .* Pol );
      P(i,k,3) = integral_dkxy( kp2 .* Pol );
      P(i,k,4) = sum(P(i,k,1:3));
    else
      P(i,k) = 2 * integral_dkxy(kern);
    end
  end
end

%% Multiply with B^2/v/L^(1-cb_in(1)) for PSD in nT^2/Hz
P = B^2 / v/L^(1-cb_in(1)) .* P;