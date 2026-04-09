import argparse
import math
import sys

# 2019 CODATA recommended value
AVOGADRO = 6.02214076e23

# nm^3 per litre  (1 L = 1 dm^3 = (10 nm)^3 = 1e24 nm^3 ... wait:
#   1 nm = 1e-9 m  →  1 nm^3 = 1e-27 m^3
#   1 L  = 1e-3 m^3  →  1 L = 1e-3 / 1e-27 nm^3 = 1e24 nm^3)
NM3_PER_LITRE = 1e24


def box_from_n_mol_and_conc(n_mol, conc):
    """Return the side (nm) of a cubic box holding *n_mol* molecules at *conc* M."""
    volume_L = n_mol / (conc * AVOGADRO)
    volume_nm3 = volume_L * NM3_PER_LITRE
    return volume_nm3 ** (1 / 3)


def conc_from_n_mol_and_volume(n_mol, volume_nm3):
    """Return the molar concentration for *n_mol* molecules in *volume_nm3* nm³."""
    volume_L = volume_nm3 / NM3_PER_LITRE
    moles = n_mol / AVOGADRO
    return moles / volume_L


def n_mol_from_conc_and_volume(conc, volume_nm3):
    """Return the number of molecules for concentration *conc* M in *volume_nm3* nm³."""
    volume_L = volume_nm3 / NM3_PER_LITRE
    return int(round(conc * AVOGADRO * volume_L))


def sphere_volume(radius_nm):
    return (4 / 3) * math.pi * radius_nm**3


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Interconvert number of molecules, molar concentration, and box geometry.\n\n"
            "Provide exactly two of the three quantities (n_mol, conc, volume/geometry)\n"
            "and the third will be computed.  Passing a geometry argument alone prints\n"
            "the corresponding volume."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--n_mol",
        type=int,
        help="Number of molecules",
    )
    parser.add_argument(
        "--conc",
        type=float,
        help="Concentration in Molar (M)",
    )
    parser.add_argument(
        "--volume",
        type=float,
        help="Volume in nm³",
    )
    parser.add_argument(
        "--sphere_r",
        type=float,
        help="Radius of a sphere in nm (converted to volume)",
    )
    parser.add_argument(
        "--cubic_side",
        type=float,
        help="Side length of a cube in nm (converted to volume)",
    )

    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        sys.exit(0)

    # --- validate positivity of anything that was supplied ---
    for name, value in vars(args).items():
        if value is not None and value <= 0:
            parser.error(f"--{name} must be a positive number (got {value})")

    # --- resolve geometry arguments to a single volume ---
    volume = args.volume
    if args.sphere_r is not None:
        if volume is not None:
            parser.error("--sphere_r and --volume are mutually exclusive")
        volume = sphere_volume(args.sphere_r)
        print(f"sphere volume: {volume:.4f} nm³  (cubic-equivalent side: {volume**(1/3):.5f} nm)")
    if args.cubic_side is not None:
        if volume is not None:
            parser.error("--cubic_side and --volume (or --sphere_r) are mutually exclusive")
        volume = args.cubic_side**3
        print(f"cubic volume:  {volume:.4f} nm³")

    n_mol = args.n_mol
    conc = args.conc

    # --- geometry alone: already printed above, nothing more to compute ---
    if n_mol is None and conc is None:
        if volume is None:
            parser.print_help()
        sys.exit(0)

    # --- compute the missing quantity from the two supplied ones ---
    if n_mol is not None and conc is not None:
        if volume is not None:
            parser.error("over-specified: provide at most two of (--n_mol, --conc, volume/geometry)")
        side = box_from_n_mol_and_conc(n_mol, conc)
        print(f"cubic box side: {side:.5f} nm")

    elif n_mol is not None and volume is not None:
        c = conc_from_n_mol_and_volume(n_mol, volume)
        print(f"concentration:  {c:.6f} M")

    elif conc is not None and volume is not None:
        n = n_mol_from_conc_and_volume(conc, volume)
        print(f"n_mol:          {n}")

    else:
        parser.error(
            "under-specified: supply exactly two of (--n_mol, --conc, " "--volume / --sphere_r / --cubic_side)"
        )
