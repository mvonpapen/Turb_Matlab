% ========================================================================
%% Title: Forward Modeling of Power Spectra (Parallelized Version)
%  Authors: Michael von Papen, Joachim Saur
%  Date: 19.03.2015
% ========================================================================
%
% This code calculates one-dimensional reduced power spectra P for given
% frequencies f and field-to-flow angles theta from a critically balanced
% energy distribution in three-dimensional k-space for MHD, ion and
% electron kinetic scales.
%
%   P = cb_psd_par (f, theta, varargin)
%
%
%% Description of input data:
%      f       Frequency in Hz
%      theta   Field-to-flow angle in degrees
%
%% Optional input data:
%      fun     function f(u) to evaluate ( 'exp', 'expdamp', 'gauss', 'gaussdamp' )
%      B       Magnetic field strength in nT
%      bounds  Boundary in log10 for integration ( bounds = [kmin kmax] )
%      cb      Critical balance exponent in [MHD KAW] range
%      L       Outer scale in m
%      n       Number of nodes in one dimension
%      ratioTP Ratio between toroidal T and poloidal fluctuations P
%      rhoi    Ion gyro radius (or other controlling kinetic ion scale) in m
%      rhoe    Electron gyro radius (or else)
%      si      Spectral index in [MHD kinetic] range
%      v       Relative plasma velocity in m/s
%      va      Alfven velocity in m/s
%
% 
%% Description of output data:
%      P       One-dimensional reduced power spectrum in nT^2/Hz
%              P is a nf x ntheta x m matrix, where nf is number of
%              frequencies, ntheta number of field-to-flow angles. If
%              ratioTP is zero, m=1 and P is the trace of the spectral
%              tensor. If ratioTP>0, m=4 and P=[Pxx Pyy Pzz Psum], where
%              Psum is the trace of the spectral tensor and Pii the
%              diagonal components.
% 
%
%% Additional comments:
% Please acknowledge the use of this code in any publications.
% Reference: von Papen, M., Saur, J., Forward Modeling of Power Spectra
%            From Three-Dimensional k-Space, ApJ, 2015
%
% More information can also be found in:
%   von Papen, M., Turbulence in Saturn's Magnetosphere and Forward
%   Modeling of Reduced Spectra from Three-Dimensional Wave Vector Space,
%   PhD thesis, Cologne, Germany, 2014
%
% I would be grateful to receive a copy of such publications under:
%  Dr. Michael von Papen
%  Institute of Geophysics
%  Pohligstr. 3
%  50969 Cologne, Germany
%  E-mail: vonpapen@geo.uni-koeln.de
%
% Updates to this code will be published under:
% www.geomet.uni-koeln.de/en/research/turbulence


function P = cb_psd_par (f, theta, varargin)


inPars = inputParser;

%% Check Input
addOptional(inPars, 'fun', 'expdamp', @(fn) any(validatestring( fn ...
                      , {'exp','expdamp','gauss','gaussdamp'})));
addOptional(inPars, 'B', 5, @isnumeric);
addOptional(inPars, 'bounds', [-10 -2], @(bn) ...
                       isequal(size(bn), [1 2]) && bn(1)<bn(2));
addOptional(inPars, 'cb', [2/3 1/3], @isnumeric);
addOptional(inPars, 'L', 1e9, @isnumeric);
addOptional(inPars, 'n', 500, @isnumeric);
addOptional(inPars, 'ratioTP', 0, @isnumeric);
addOptional(inPars, 'rhoi', 1e5, @isnumeric);
addOptional(inPars, 'rhoe', NaN, @isnumeric);
addOptional(inPars, 'si', [-10/3 -11/3], @isnumeric);
addOptional(inPars, 'v', 6e5, @isnumeric);
addOptional(inPars, 'va', 6e4, @isnumeric);

parse(inPars, varargin{:});

fun = inPars.Results.fun;
L = inPars.Results.L;
bounds = inPars.Results.bounds;
cb = inPars.Results.cb;
B = inPars.Results.B;
rhoi = inPars.Results.rhoi;
if isnan(inPars.Results.rhoe);
       rhoe = rhoi/42.85;
  else rhoe = inPars.Results.rhoe;
end
va = inPars.Results.va;
v = inPars.Results.v;
ratioTP = inPars.Results.ratioTP;
n = inPars.Results.n;
si = inPars.Results.si;
ntheta=length(theta);


%% Wave vector in rotated coordinate system
kix = 2*pi*f/v;
nf=length(kix);

%% Basic Parameters
si1 = si(1); %-10/3 -> k^{-5/3}
si2 = si(2); %-11/3 -> k^{-7/3}
cb1 = cb(1); % 2/3  -> Alfven
cb2 = cb(2); % 1/3  -> KAW


%% k-space gridpoints
kmin=bounds(1);
kmax=bounds(2);

% Range to cover in each dimension of the k-integration.
% n is the nunber of nodes integrated over for each sign in each direction,
% kiz (=kz') covers both signs, there are thus 2n nodes in that direction.
kRange = 10.^(kmin + (0:n-1)*(kmax-kmin)/n);

% Log-distribute z-nodes over both positive and negative range.
kiz = repmat([-kRange(end:-1:1) kRange], n, 1);

% The y-component only needs positive values, due to symmetry.
% ky' is the same as ky, hence we don't use the name kiy here.
ky = repmat(transpose(kRange), 1, 2*n);



%% Ion-cyclotron frequency
% Ion-cyclotron damping sets in at k_para ~ w_ic/V_A
% Thus, parallel scales cannot reach k_para >> w_ic/V_A
wic = 1.6e-19*B*1e-9/1.67e-27; %wic for protons


%% Set output variable
if ratioTP ~=0
    Pxx = zeros(nf,ntheta);
    Pyy = zeros(nf,ntheta);
    Pzz = zeros(nf,ntheta);
else
    P = zeros(nf,ntheta);
end


%% Begin with loop over theta
parfor k=1:ntheta

  if ratioTP ~=0
      Pxxtmp = zeros(nf,1);
      Pyytmp = zeros(nf,1);
      Pzztmp = zeros(nf,1);
  else
      Ptmp = zeros(nf,1);
  end
  
  
  %% Begin loop over frequency
  for i=1:nf

      
    %% Calculate wave numbers in unprimed coordinate system (see Eqs 6-8)
    %% for the sake of brevity
    
    % CAUTION: Problems may arise when setting kiz=0 on x-axis for small
    % outer scales and/or small angles on MHD scales. In that case comment
    % lines 167 & 171-177
    if theta(k)==0
      kx = kix(i)*sind(theta(k)) - kiz.*cosd(theta(k)); %=kx in unrotated system
      %ky=kiy in unrotated system
      kz = kix(i)*cosd(theta(k)) + kiz.*sind(theta(k)); %=kz in unrotated system
    else
      kx = kix(i)*sind(theta(k)) - kiz.*cosd(theta(k))                             ...
          + kix(i).*cosd(theta(k))^2/sind(theta(k)); %=kx in unrotated system with kiz=0 on x-axis
      %ky=kiy in unrotated system
      kz = kiz.*sind(theta(k)); %=kz in unrotated system with kiz=0 on x-axis
    end

    kern = zeros(n, 2*n);

    kp2 = ky.^2 + kx.^2; % k_perp^2
    kabs2 = kp2 + kz.^2; % |k|^2

    
    
    %% Single components Alfven cascade.
    
    i1 = find(kp2 <= 1/rhoi^2 & kabs2 > 1/L^2);
    
    switch fun
        
      case 'exp' % Eq. (14)
          kern(i1) = kp2(i1).^(si1/2)                                                       ...
                         .* exp(-L^(1-cb1) .* abs( kz(i1) ) ./ kp2(i1).^(cb1/2) );
      case 'expdamp' % Eqs. (14), (20), (22)
          kern(i1) = kp2(i1).^(si1/2)                                                       ...
                         .* exp( - L^(1-cb1) .* abs(kz(i1) ./ kp2(i1).^(cb1/2))              ...
                                 - sqrt( kp2(i1) ) * rhoe                                  ...
                                 - abs(kz(i1)) * va / wic );
      case 'gauss' % f(u) given by Gaussian centered at critical balance
          kern(i1) = kp2(i1).^(si1/2)                                                       ...
                         .* exp( - (L^(1-cb1) * abs( kz(i1) ) ./ kp2(i1).^(cb/2) - 1 ).^2 ) ...
                         ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = kp2(i1).^(si1/2)                                                       ...
                         .* exp( - ( L^(1-cb1) * abs(kz(i1)) ./ kp2(i1).^(cb1/2) - 1 ).^2    ...
                                 - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )         ...
                         ./ sqrt(pi); 
                     
    end

    

    %% Single components KAW cascade
    
    i1 = find(kp2 > 1/rhoi^2 & kabs2 > 1/L^2 & kp2<= 1/rhoe^2);
    
    switch fun
        
      case 'exp' % Eq. (16)
          kern(i1) = rhoi^(si2-si1)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( -L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) );
      case 'expdamp' % Eqs. (16), (20), (22)
          kern(i1) = rhoi^(si2-si1)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( - L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2)            ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic );
      case 'gauss' % f(u) given by Gaussian centered at critical balance
          kern(i1) = rhoi^(si2-si1)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( -(L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) - 1 ).^2 ) ...
                        ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = rhoi^(si2-si1)                                                                      ...
                        .* kp2(i1).^(si2/2)                                                            ...
                        .* exp( - (L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) ./ kp2(i1).^(cb2/2) - 1 ).^2  ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )                      ...
                        ./ sqrt(pi);
                    
    end

    
    
    %% Single components of cascade at electron scales
    
    i1 = find(kp2 > 1/rhoe^2 & kabs2 > 1/L^2);
    
    switch fun
        
      case 'exp' % Eq. (18)
          kern(i1) = rhoi^(si2-si1)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( -L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) * rhoe^cb2 );
      case 'expdamp' % Eqs. (18), (20), (22)
          kern(i1) = rhoi^(si2-si1)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( - L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) * rhoe^cb2            ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic);
      case 'gauss' % f(u) given by Gaussian centered at critical balance
          kern(i1) = rhoi^(si2-si1)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( -(L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) * rhoe^cb2 - 1 ).^2 ) ...
                        ./ sqrt(pi);
      case 'gaussdamp'
          kern(i1) = rhoi^(si2-si1)                                                             ...
                        .* kp2(i1).^(si2/2)                                                   ...
                        .* exp( - (L^(1-cb1) * rhoi^(cb1-cb2) * abs(kz(i1)) * rhoe^cb2 - 1 ).^2  ...
                                - sqrt(kp2(i1)) * rhoe - abs(kz(i1)) * va / wic )             ...
                        ./ sqrt(pi);
                    
    end

    
    
    %% Algorithm for 2D-integration over the kxy-domain ranges.
    integral_dkxy = @(fk) trapz( ky(:,1), trapz(kiz(1,:), fk, 2) );

    
    
    %% Calculate the components/trave of spectral tensor
    
    if ratioTP ~= 0
        
      Tor = ratioTP / (1+ratioTP) * kern ./ kp2;      
      Pol = kern / ratioTP ./ kabs2;

      % Diagonal components of spectral tensor, Eqs. (1-3)
      Pxxtmp(i) = integral_dkxy( ky.^2.*Tor + (kx.*kz).^2 ./ kp2 .* Pol );
      Pyytmp(i) = integral_dkxy( kx.^2.*Tor + (ky.*kz).^2 ./ kp2 .* Pol );
      Pzztmp(i) = integral_dkxy( kp2 .* Pol );
      
    else
    
      Ptmp(i) = integral_dkxy(kern);
      
    end
    
  end
  
  if ratioTP==0
      P(:,k)=Ptmp;
  else
      Pxx(:,k)=Pxxtmp;
      Pyy(:,k)=Pxxtmp;
      Pzz(:,k)=Pxxtmp;
  end
  
end

if ratioTP~=0
    P(:,:,4)=Pxx+Pyy+Pzz;
    P(:,:,3)=Pzz;
    P(:,:,2)=Pyy;
    P(:,:,1)=Pxx;
end

%% Multiply with B^2/v/L^(1-cb1) for PSD in nT^2/Hz (Eq. (11))
% also multiply with 2 because of symmetry in ky
P = 2 * B^2 / v / L^(1-cb1) .* P;