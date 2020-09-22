function I = testKramersKronig(omega,chi)
  % Test whether data fulfills the Kramers-Kronig relation.
  % 
  % I = testKramersKronig(omega,chi)
  %   omega is a vector of frequencies
  %   chi is a vector of complex susceptibility data
  %   I is a vector of the real part of the susceptibility, determined from
  %   the imaginary part of chi
  % 
  % I is calculated by numerically integrating the imaginary part of chi.
  % Is is then plotted alongside real(chi) and -imag(chi) to compare.
  % 
  % The integration is done using Maclaurin's formula from as described in:
  % Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
  % Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).

  I = integrateMaclaurin(omega,-imag(chi)); % Do the calculation

  figure(91); clf; hold on; box on; grid on; % Plot on a figure that's not in use
  plot( ...
    log10(omega), ...
    real(chi), ...
    'LineWidth',2,'DisplayName','χ′')
  plot( ...
    log10(omega), ...
    -imag(chi), ...
    'LineWidth',2,'DisplayName','χ′′')
  plot( ...
    log10(omega), ...
    I, ...
    'LineWidth',2,'DisplayName','I')
  legend
  xlabel('ω')
  ylabel('χ')
  setPlotSize(15,10,'centimeters')
end