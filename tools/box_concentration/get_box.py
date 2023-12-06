import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO!")
    parser.add_argument("--n_mol", type=int, required=True, help="Number of molecules to simulate.")
    parser.add_argument("--conc", type=float, required=True, help="Concentration in Molar")

    args = parser.parse_args()

    n_mol = args.n_mol
    conc = args.conc
    box = 10**8 * (n_mol / (conc * 6.022 * 10**23)) ** (1 / 3)
    print("cubic box of side: %.5f nm" % box)
