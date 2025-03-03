import numpy as np
from vaspwfc import vaspwfc
import sys

def read_wavecar(filename='WAVECAR'):
    """Read eigenvalues and basic info from WAVECAR."""
    try:
        wav = vaspwfc(filename)
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        sys.exit(1)
    
    nkpts = wav._nkpts  # Number of k-points
    nbands = wav._nbands  # Total number of bands
    nspin = wav._nspin  # Number of spins
    efermi = wav._efermi  # Fermi level
    
    # Read eigenvalues (shape: nspin, nkpts, nbands)
    energies = wav._bands
    occupations = wav._occs 
    
    return energies, occupations, efermi, nkpts, nbands, nspin

def degeneracy_check(filenames, tol_degeneracy=1e-4):
    """Check degeneracy-allowed numbers of bands from WAVECAR files."""
    minband = None
    minvband = None
    mincband = None
    ok = None
    ok_vband = None
    ok_cband = None
    
    for ifile, infile in enumerate(filenames):
        # Read data from WAVECAR
        energies, occupations, efermi, nkpts, nbands, nspin = read_wavecar(infile)
        
        # Print basic info
        print(f"\nReading eigenvalues from file {infile}")
        print(f"Number of spins:     {nspin}")
        print(f"Number of bands:     {nbands}")
        print(f"Number of k-points:  {nkpts}")
        
        # Determine valence and conduction bands based on Fermi level
        # Assume fully occupied (occ > 0.5) or below efermi for valence
        ifmax = np.zeros((nspin, nkpts), dtype=int)  # Index of highest occupied band
        for ispin in range(nspin):
            for ik in range(nkpts):
                # Find highest band below efermi or with significant occupation
                occ = occupations[ispin, ik, :]
                for ib in range(nbands):
                    if energies[ispin, ik, ib] >= efermi and occupations[ispin, ik, ib] < 0.5:
                       ifmax[ispin, ik] = ib - 1
                       continue
        
        nvband = min(np.max(ifmax, axis=1) + 1)  # Number of valence bands (counting from 0)
        ncband = nbands - max(np.max(ifmax, axis=1) + 1)  # Number of conduction bands
        
        # Initialize for first file
        if ifile == 0:
            minband = nbands
            minvband = nvband
            mincband = ncband
            
            ok = np.ones(minband - 1, dtype=bool)
            ok_vband = np.ones(minvband - 1, dtype=bool)
            if mincband - 1 > 0:
                ok_cband = np.ones(mincband - 1, dtype=bool)
            else:
                ok_cband = np.array([], dtype=bool)
        else:
            minband = min(minband, nbands)
            minvband = min(minvband, nvband)
            mincband = min(mincband, ncband)
        
        # Check total bands
        for ib in range(minband - 1):
            energy_diff = np.abs(energies[:, :, ib] - energies[:, :, ib + 1])
            ok[ib] = ok[ib] & np.all(energy_diff > tol_degeneracy)
        
        # Check valence bands (counting down from ifmax)
        for ib in range(minvband - 1):
            for ispin in range(nspin):
                for ik in range(nkpts):
                    idx_high = ifmax[ispin, ik] - ib
                    idx_low = ifmax[ispin, ik] - ib - 1
                    energy_diff = abs(energies[ispin, ik, idx_high] - energies[ispin, ik, idx_low])
                    ok_vband[ib] = ok_vband[ib] & (energy_diff > tol_degeneracy)
        
        # Check conduction bands (counting up from ifmax)
        for ib in range(mincband - 1):
            for ispin in range(nspin):
                for ik in range(nkpts):
                    idx_low = ifmax[ispin, ik] + ib + 1
                    idx_high = ifmax[ispin, ik] + ib + 2
                    energy_diff = abs(energies[ispin, ik, idx_low] - energies[ispin, ik, idx_high])
                    ok_cband[ib] = ok_cband[ib] & (energy_diff > tol_degeneracy)
    
    # Output results
    if len(filenames) > 1:
        print(f"\nMinimum number of bands in files: {minband}")
    
    print("\n== Degeneracy-allowed numbers of bands (for epsilon and sigma) ==")
    for ib in range(minband - 1):
        if ok[ib]:
            print(ib + 1)  # +1 to match 1-based indexing in output
    print(f"Note: cannot assess whether or not highest band {minband} is degenerate.")
    
    print("\n== Degeneracy-allowed numbers of valence bands (for inteqp, kernel, and absorption) ==")
    for ib in range(minvband - 1):
        if ok_vband[ib]:
            print(ib + 1)
    print(minvband)  # Using all valence bands is always allowed
    
    if mincband - 1 > 0:
        print("\n== Degeneracy-allowed numbers of conduction bands (for inteqp, kernel, and absorption) ==")
        for ib in range(mincband - 1):
            if ok_cband[ib]:
                print(ib + 1)
        print(f"Note: cannot assess whether or not highest conduction band {mincband} is degenerate.")

def main():
    # Get command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python degeneracy_check.py WAVECAR [WAVECAR2 ...]")
        sys.exit(1)
    
    filenames = sys.argv[1:]
    degeneracy_check(filenames, tol_degeneracy=1e-4)

if __name__ == "__main__":
    main()
