import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate number of moles, concentration, box size or volume based on provided parameters."
    )
    parser.add_argument(
        "--n_mol",
        type=int,
        required=False,
        help="Number of molecules, to be combined with --conc or one from (--volume, --sphere_r, --cubic_side)",
    )
    parser.add_argument(
        "--conc",
        type=float,
        required=False,
        help="Concentration in Molar, to be combined with --n_mol or one from (--volume, --sphere_r, --cubic_side)",
    )
    parser.add_argument("--volume", type=float, required=False, help="Volume in nm^3, to be combined with --n_mol or --conc")
    parser.add_argument("--sphere_r", type=float, required=False, help="Radius of a sphere in nm")
    parser.add_argument("--cubic_side", type=float, required=False, help="Side of a cube in nm")

    args = parser.parse_args()

    # Print help if no arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        exit(0)

    n_mol = args.n_mol
    conc = args.conc
    volume = args.volume
    sphere_r = args.sphere_r
    cubic_side = args.cubic_side

    if n_mol is not None and n_mol > 0 and conc is not None and conc > 0:
        box = 10**8 * (n_mol / (conc * 6.022 * 10**23)) ** (1 / 3)
        print("cubic box of side: %.5f nm" % box)
        exit()

    if sphere_r is not None and sphere_r > 0:
        volume = (4 / 3) * np.pi * sphere_r**3

    if cubic_side is not None and cubic_side > 0:
        volume = cubic_side**3

    if n_mol is not None and n_mol > 0 and volume is not None and volume > 0:
        conc = 10**24 * (n_mol / (volume * 6.022 * 10**23))
        print("concentration is: %.12f M" % conc)
        exit()

    if conc is not None and conc > 0 and volume is not None and volume > 0:
        n_mol = int(round((volume * conc * 6.022 * 10**23) / 10**24))
        print("n_mol is: %i" % n_mol)
        exit()

    if sphere_r is not None and sphere_r > 0:
        volume = (4 / 3) * np.pi * sphere_r**3
        side = volume ** (1 / 3)
        print("volume is %.6f nm^3 and cubic side would be %.6f nm" % (volume, side))
        exit()

    parser.print_help()
