"""
ShearFunctions: Python translation of MATLAB functions for shear data analysis.

Functions for numerical integration of imaginary parts of complex susceptibility
data using the Kramers-Kronig relation and Maclaurin's formula.

Reference:
Ohta, K. & Ishida, H. Comparison among Several Numerical Integration Methods 
for Kramers-Kronig Transformation. Applied Spectroscopy (1988).
"""

import numpy as np
import matplotlib.pyplot as plt


def integrate_maclaurin(omega, chi_imag):
    """
    Use Maclaurin's formula to integrate the imaginary part of complex 
    susceptibility data to obtain the real part.
    
    Parameters
    ----------
    omega : array-like
        Frequency array
    chi_imag : array-like
        Imaginary part of complex susceptibility data
    
    Returns
    -------
    I : ndarray
        Integrated real part of susceptibility
    
    Notes
    -----
    The integration is done using Maclaurin's formula as described in:
    Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
    Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).
    """
    omega = np.asarray(omega)
    chi_imag = np.asarray(chi_imag)
    
    N = len(chi_imag)
    d_omega = np.diff(omega)
    
    # Set up variables for the two cases in the loop below
    indices_even = np.arange(0, N, 2)  # 0-indexed: 0, 2, 4, ...
    indices_odd = np.arange(1, N, 2)   # 0-indexed: 1, 3, 5, ...
    d_omega_even = d_omega[0::2]
    d_omega_odd = d_omega[1::2]
    
    # One of the vectors will be one value too short, so we lengthen it 
    # by repeating the last value according to a formula
    if N % 2 == 0:  # Even number of points
        d_omega_odd = np.append(d_omega_odd, d_omega_odd[-1]**2 / d_omega_odd[-2])
    else:  # Odd number of points
        d_omega_even = np.append(d_omega_even, d_omega_even[-1]**2 / d_omega_even[-2])
    
    I = np.full(chi_imag.shape, np.nan)
    
    for index_current in range(N):
        if index_current % 2 == 1:  # Even index in MATLAB (1-indexed)
            f = (omega[indices_even] * chi_imag[indices_even]) / \
                (omega[indices_even]**2 - omega[index_current]**2)
            hf = f * d_omega_even
        else:  # Odd index in MATLAB (1-indexed)
            f = (omega[indices_odd] * chi_imag[indices_odd]) / \
                (omega[indices_odd]**2 - omega[index_current]**2)
            hf = f * d_omega_odd
        
        I[index_current] = np.sum(hf)
    
    I = I * 4 / np.pi
    
    return I


def integrate_g_kramers_kronig(s, shear_name, fres):
    """
    Numerically integrate the imaginary part of inverted shear data using 
    the Kramers-Kronig relation.
    
    Parameters
    ----------
    s : dict
        Dictionary containing inverted shear data with structure:
        s[shear_name]['fr'] : array
            Frequency array (can be 2D or 3D)
        s[shear_name]['G'] : array
            Complex shear modulus (same dimensions as fr)
    shear_name : str
        Name of the inversion method whose results to integrate
    fres : float
        Approximate resonance frequency of the system
        (frequencies above this are excluded from integration)
    
    Returns
    -------
    s : dict
        Updated dictionary with two new fields added to s[shear_name]:
        - 'KramersKronig': array of integration results
        - 'KramersKronigShift': mean of first few values (for offset correction)
    
    Notes
    -----
    Comparing s[shear_name]['KramersKronig'] to real(s[shear_name]['G']) 
    gives an indication of how well data fulfills the Kramers-Kronig relation.
    
    The integration is done using Maclaurin's formula as described in:
    Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
    Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).
    """
    fr = np.asarray(s[shear_name]['fr'])
    G = np.asarray(s[shear_name]['G'])
    
    # Find indices where frequency is below resonance
    if fr.ndim == 3:
        fr_end = fr[:, -1, -1]
    elif fr.ndim == 2:
        fr_end = fr[:, -1]
    else:
        fr_end = fr
    
    fr_indices = np.where(fr_end < fres)[0]
    
    # Get dimensions (handles both 2D and 3D arrays)
    shape = fr.shape
    m1, m2 = shape[0], shape[1]
    m3 = shape[2] if len(shape) == 3 else 1
    
    # Initialize output arrays
    kramers_kronig = np.full(shape, np.nan)
    kramers_kronig_shift = np.full((m3, m2), np.nan)
    
    # Perform integration for each layer and column
    for index_layer in range(m3):
        for index_column in range(m2):
            if fr.ndim == 3:
                fr_data = fr[fr_indices, index_column, index_layer]
                g_data = G[fr_indices, index_column, index_layer]
            else:
                fr_data = fr[fr_indices, index_column]
                g_data = G[fr_indices, index_column]
            
            I = integrate_maclaurin(fr_data, -np.imag(g_data))
            
            if fr.ndim == 3:
                kramers_kronig[fr_indices, index_column, index_layer] = I
            else:
                kramers_kronig[fr_indices, index_column] = I
            
            # Calculate shift as mean of first 5 negative values
            negative_indices = np.where(I < 0)[0][:5]
            if len(negative_indices) > 0:
                kramers_kronig_shift[index_layer, index_column] = np.mean(I[negative_indices])
    
    s[shear_name]['KramersKronig'] = kramers_kronig
    s[shear_name]['KramersKronigShift'] = kramers_kronig_shift
    
    return s


def test_kramers_kronig(omega, chi, show_plot=True):
    """
    Test whether data fulfills the Kramers-Kronig relation.
    
    Parameters
    ----------
    omega : array-like
        Vector of frequencies
    chi : array-like
        Vector of complex susceptibility data
    show_plot : bool, optional
        If True, display a plot comparing the real part, imaginary part,
        and integrated values. Default is True.
    
    Returns
    -------
    I : ndarray
        Vector of the real part of susceptibility, determined from
        the imaginary part of chi via Kramers-Kronig transformation
    
    Notes
    -----
    The integration is done using Maclaurin's formula as described in:
    Ohta, K. & Ishida, H. Comparison among Several Numerical Integration
    Methods for Kramers-Kronig Transformation. Applied Spectroscopy (1988).
    """
    omega = np.asarray(omega)
    chi = np.asarray(chi)
    
    # Do the calculation
    I = integrate_maclaurin(omega, -np.imag(chi))
    
    if show_plot:
        plt.figure(91, figsize=(10, 6))
        plt.clf()
        plt.grid(True, alpha=0.3)
        
        log_omega = np.log10(omega)
        
        plt.plot(log_omega, np.real(chi), linewidth=2, label="χ' (real part)")
        plt.plot(log_omega, -np.imag(chi), linewidth=2, label="χ'' (-imaginary part)")
        plt.plot(log_omega, I, linewidth=2, label="I (integrated)")
        
        plt.legend(fontsize=11)
        plt.xlabel('log₁₀(ω)', fontsize=12)
        plt.ylabel('χ', fontsize=12)
        plt.title('Kramers-Kronig Test: Real part vs Integration Result', fontsize=13)
        plt.tight_layout()
        plt.show()
    
    return I


# Example usage and testing
if __name__ == "__main__":
    print("=" * 70)
    print("ShearFunctions: Kramers-Kronig Transformation Test")
    print("=" * 70)
    
    # Create synthetic test data
    # Use a Lorentzian-like response for testing
    omega = np.logspace(-2, 2, 100)  # Frequencies from 0.01 to 100
    
    # Synthetic complex susceptibility data (Lorentzian model)
    omega_0 = 1.0  # Resonance frequency
    gamma = 0.1    # Damping
    chi_real_test = 1.0 / ((omega_0**2 - omega**2) / omega_0**2 + 1j * gamma * omega / omega_0)
    
    print(f"\nTest data created:")
    print(f"  Frequency range: {omega[0]:.4f} to {omega[-1]:.4f}")
    print(f"  Number of points: {len(omega)}")
    print(f"  Resonance frequency: {omega_0}")
    print(f"  Damping: {gamma}")
    
    # Run the test
    print(f"\nRunning Kramers-Kronig test...")
    I = test_kramers_kronig(omega, chi_real_test, show_plot=True)
    
    # Compare results
    chi_real_actual = np.real(chi_real_test)
    error = np.abs(chi_real_actual - I)
    relative_error = error / (np.abs(chi_real_actual) + 1e-10)
    
    print(f"\nResults:")
    print(f"  Mean absolute error: {np.nanmean(error):.6e}")
    print(f"  Max absolute error: {np.nanmax(error):.6e}")
    print(f"  Mean relative error: {np.nanmean(relative_error):.6e}")
    print(f"  Max relative error: {np.nanmax(relative_error):.6e}")
    print(f"\nThe Kramers-Kronig relation is satisfied when I ≈ χ' (real part)")
