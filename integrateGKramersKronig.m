function s = integrateGKramersKronig(s,shearName,fres)
  % Function to numerically integrate the imaginary part of inverted shear
  % data using the Kramers-Kronig relation.
  % 
  % s = integrateGKramersKronig(s,shearName,fres)
  %   s is a structure of inverted shear data
  %   shearName is a string contining the name of the inversion methods
  %     whose results you want to integrate
  %   fres is an approximate value for the resonance frequency of the PSG
  % 
  % The function adds two fields to the structure s:
  %   s.(shearName).KramersKronig is an array of integration results of the
  %     same dimensions as s.(shearName).fr
  %   s.(shearName).KramersKronigShift is an array containing the mean of
  %     the first few values of s.(shearName).KramersKronig at each
  %     temperature. This may be used as a substitute for the integration
  %     constant when plotting the results
  % 
  % Comparing s.(shearName).KramersKronig to real(s.(shearName).G) gives an
  % indication of how well data fulfills the Kramers-Kronig relation.
  % 
  % The integration is done using Maclaurin's formula from as described in:
  % Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
  % Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).

  frIndices = find(s.(shearName).fr(:,end,end) < fres); % If the resonance is included, the integral is completely wrong
  [m1,m2,m3] = size(s.(shearName).fr); % This ensures that the functions accepts G-arrays of both two and three dimensions
  
  s.(shearName).KramersKronig = NaN(m1,m2,m3);
  s.(shearName).KramersKronigShift = NaN(m3,m2);

  for indexLayer = 1:m3
    for indexColumn = 1:m2
      I = integrateMaclaurin(s.(shearName).fr(frIndices,indexColumn,indexLayer),-imag(s.(shearName).G(frIndices,indexColumn,indexLayer)));
      s.(shearName).KramersKronig(frIndices,indexColumn,indexLayer) = I;
      s.(shearName).KramersKronigShift(indexLayer,indexColumn) = mean(I(find(I < 0,5,'first')));
    end
  end
end