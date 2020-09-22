function I = integrateMaclaurin(omega,chiImag)
  % Use Maclaurin's formula to integrate the imaginary part of complex 
  % susceptibility data to obtain the real part.
  % 
  % I = integrateMaclaurin(omega,chiImag)
  % 
  % The integration is done using Maclaurin's formula from as described in:
  % Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
  % Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).
  
  N = numel(chiImag);
  dOmega = diff(omega);

  % Set uo variables for the two cases in the loop below
  indicesEven = 1:2:N;
  indicesOdd = 2:2:N;
  dOmegaEven = dOmega(1:2:end);
  dOmegaOdd = dOmega(2:2:end);

  if mod(N,2) == 0 % One of the vectors will be one value too short, so we lengthen it by repeating the last value:
    dOmegaOdd(end + 1) = dOmegaOdd(end)^2/dOmegaOdd(end - 1);
  else
    dOmegaEven(end + 1) = dOmegaEven(end)^2/dOmegaEven(end - 1);
  end

  I = NaN(size(chiImag));
  
  for indexCurrent = 1:N
    if mod(indexCurrent,2) == 0
      f = (omega(indicesEven).*chiImag(indicesEven))./(omega(indicesEven).^2 - omega(indexCurrent)^2);
      hf = f.*dOmegaEven;
    else
      f = (omega(indicesOdd).*chiImag(indicesOdd))./(omega(indicesOdd).^2 - omega(indexCurrent)^2);
      hf = f.*dOmegaOdd;
    end

    I(indexCurrent) = sum(hf);
  end

  I = I*4/pi;
end