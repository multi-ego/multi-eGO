from fileinput import filename

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
from matplotlib.cm import ScalarMappable

bkbn_attypes = ["H", "O", "N", "C", "CAH"]#, "CAH2"]
H_interactions = ["O", "H", "OA", "OM", "OS"]

noplot = ["C0", "CH3p", "CR1", "OE", "P"]

class AttrDict(dict):
    def __getattr__(self, key):
        return self[key]

# TODO re-define the pickle to avoid the utils module, which is not available anymore. This is a temporary fix to load the pickle file.
class RemapUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "utils":
            # point this at wherever that class/function now lives
            module = __name__  # this file's current module name
        return super().find_class(module, name)

class InteractionMatrix:

    def __init__(self, pkl_file, emax=0.16, c12_rep_df=None, pth=None, show=False, pth_bkbn=None, emax_bkbn=None, pth_sc_bkbn=None, emax_sc_bkbn=None, om_scaledown=None, pkl_file2=None, r2_corr_sb=None):
        self.atmat = self.build_matrix_dataframe(pkl_file, filename2=pkl_file2, r2_corr_sb=r2_corr_sb)
        self.emax = emax
        self.emax_bkbn = emax_bkbn if emax_bkbn is not None else emax
        self.emax_sc_bkbn = emax_sc_bkbn if emax_sc_bkbn is not None else self.emax_bkbn
        self.pth = pth
        self.pth_bkbn = pth_bkbn if pth_bkbn is not None else pth
        self.pth_sc_bkbn = pth_sc_bkbn if pth_sc_bkbn is not None else (self.pth + self.pth_bkbn) / 2
        self.om_scaledown = om_scaledown
        self.r2_corr_sb = r2_corr_sb
        self.c12_rep_df = c12_rep_df

        # self.define_p_threshold()
        self.nonlocal_matrix()
        if show:
            self.plot_energy_matrix(show=show)

    def roundPartial(self, value, resolution=0.005):
        """
        Regularize the epsilon values to have minimum resolution in decimal places
        """
        return float(round (value / resolution) * resolution)

    def regall(self, ps, preference, emax, pth):
        """
        Regularize the energy values to have minimum resolution in decimal places
        Arguments:
        ps: probability of the atomtype pair
        preference: (Now is 1). In principle it could be a reference probability per atomtype pair to reweight probabilities (e.g. water interaction probabilities)
        emax: maximum energy value
        pth: threshold probability to split attractive and repulsive interactions
        """
        e = - emax/np.log(1e-4) * np.log((ps + 1e-10)/(np.maximum(1e-4, preference*pth)))
        e /= np.max(e)
        e = np.where(e==-np.inf, -2, e)
        e = np.where(e > 0.0, 0.07 + e* (emax-0.07), e)
        # e = np.array([self.roundPartial(ei) for ei in e])
        return e

    def plot_energy_matrix(self, out_path = None, show = False):
        """
        Plot the energy matrix as a heatmap.
        """

        # TODO read them by gromos to avoid hardcoding them here
        bkbnd_atoms = ["C", "O", "N", "CAH", "CAH2", "H"]
        attype_ordering_charged = [ "OM",  "NL", "NZ"]
        attype_ordering_polar = [ "NE", "NR", "NT", "OA","OS", "CZ", "CS", "CH1t", "CH2r", "CH2t", "CH3t", "S", "SM", "SH"]
        attype_ordering_apolar = ["CH1", "CR", "CH",  "CH2", "CH3"]



        list_atom1 = np.sort(np.array(self.atmat["attype1"].unique()))
        list_atom2 = np.sort(np.array(self.atmat["attype2"].unique()))

        # TODO read them by gromos to avoid hardcoding them here
        # remove from list_atom1 and list_atom2 the atoms in noplot
        list_atom1 = [at for at in list_atom1 if at not in noplot]
        list_atom2 = [at for at in list_atom2 if at not in noplot]

        # create ordered list of atom types based on attype_ordering
        order_list_atom1 = [at for at in bkbnd_atoms if at in list_atom1] + [at for at in attype_ordering_charged if at in list_atom1] + [at for at in attype_ordering_polar if at in list_atom1] + [at for at in attype_ordering_apolar if at in list_atom1]
        order_list_atom2 = [at for at in bkbnd_atoms if at in list_atom2] + [at for at in attype_ordering_charged if at in list_atom2] + [at for at in attype_ordering_polar if at in list_atom2] + [at for at in attype_ordering_apolar if at in list_atom2]

        energies = []
        for at1 in order_list_atom1:
            for at2 in order_list_atom2:
                if at1 in noplot or at2 in noplot:
                    continue
                energy = self.atmat.loc[((self.atmat["attype1"]==at1) & (self.atmat["attype2"]==at2)) | ((self.atmat["attype1"]==at2) & (self.atmat["attype2"]==at1)), "energy"].values
                if (at1 == "H" and at2 not in H_interactions) or (at2 == "H" and at1 not in H_interactions):
                    energies.append(-np.inf) # force repulsion for H interactions with non-interaction-H atoms
                    continue
                if len(energy) > 0:
                    energy = energy[0]
                else:
                    energy = -1
                if energy > 0:
                    energies.append(energy)
                else:
                    energies.append(-1)
                
        colors = ["tab:orange", "white", "blue"]
        cmap = LinearSegmentedColormap.from_list(
            "red_white_blue", ["tab:orange", "white", "blue"]
        )

        # Full normalization for the image
        energy_matrix = np.array(energies).reshape(len(order_list_atom1), len(order_list_atom2))
        norm = TwoSlopeNorm(vmin=-1,vcenter=0.07,   vmax=np.maximum(self.emax, self.emax_bkbn))
        # set -2 values to grey 
        energy_matrix = np.where(energy_matrix == -2, np.nan, energy_matrix)
        # set nan values to grey
        cmap.set_bad(color="grey", alpha=0.6)

        fig, ax = plt.subplots()
        im = ax.imshow(energy_matrix, cmap=cmap, norm=norm, origin="lower")
        blue_part = cmap(np.linspace(0.5, 1, 256))
        blue_cmap = LinearSegmentedColormap.from_list("blue_part", blue_part)

        pos_norm = Normalize(vmin=0.07, vmax=np.maximum(self.emax, self.emax_bkbn))
        sm = ScalarMappable(norm=pos_norm, cmap=blue_cmap)
        sm.set_array([])  # required for colorbar
        ax.set_xticks(ticks=np.arange(len(order_list_atom2)), labels=order_list_atom2, rotation=90)
        ax.set_yticks(ticks=np.arange(len(order_list_atom1)), labels=order_list_atom1)
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Attractive interaction strenght")
        ax.grid()
        if out_path is not None:
            plt.savefig(out_path)
        if show:
            plt.show()

    def read_pickle_file(self, filename):
        _this_dir = os.path.dirname(os.path.abspath(__file__))
        pikle_path = os.path.join(_this_dir, filename)

        # Fallback: older layout had the CSV in src/multiego relative to project root
        if not os.path.exists(pikle_path):
            pikle_path = os.path.join(_this_dir, "..", "src", "multiego", filename)

        with open(pikle_path, "rb") as f:
            data = RemapUnpickler(f).load()
            
        atmat = pd.DataFrame()
        atmat["atom_pair"] = data.keys()
        probs = []
        dists = []
        cutoffs = []
        bins = []
        kde = []
        tot_reps = []
        sum_probs = []
        for key in data.keys():
            # probs.append(data[key].p_repeats_n2)
            probs.append(data[key].p_density)
            dists.append(data[key].exp_aver)
            cutoffs.append(data[key].cutoff)
            bins.append(data[key].xbins)
            kde.append(data[key].kde)
            tot_reps.append(data[key].tot_repeats)
            sum_probs.append(np.sum(data[key].sum_probs))
        atmat["probability"] = np.array(probs)
        atmat["exp_aver"] = np.array(dists)
        atmat["attype1"] = atmat["atom_pair"].str.split("_").str[0]
        atmat["attype2"] = atmat["atom_pair"].str.split("_").str[1]
        atmat["cutoff"] = np.array(cutoffs)
        atmat["bins"] = bins
        atmat["kde"] = kde
        atmat["tot_repeats"] = tot_reps
        atmat["sum_probs"] = sum_probs
        # fill nan values in probability with 0
        atmat["probability"] = atmat["probability"].fillna(0)
        return atmat

    def build_matrix_dataframe(self, filename, filename2=None, r2_corr_sb=None):
        """
        Read atdhito data from pikle file per atomtype pair and return a pandas dataframe with the following columns:
        - atom_pair: atomtype pair
        - probability: probability of the atomtype pair
        - exp_aver: expected average distance of the atomtype pair
        - attype1: first atomtype
        - attype2: second atomtype
        - cutoff: cutoff distance of the atomtype pair
        """
        _this_dir = os.path.dirname(os.path.abspath(__file__))
        
        atmat = self.read_pickle_file(filename)
        
        if filename2 != None:
            atmat_sb = self.read_pickle_file(filename2)

        # use exp_aver and probability from atmat_sb for atom pairs where one of the atomtypes is in bkbn_attypes and the other is not
        if filename2 != None:
            for index, row in atmat.iterrows():
                at1 = row["attype1"]
                at2 = row["attype2"]
                if (at1 in bkbn_attypes and at2 not in bkbn_attypes) or (at2 in bkbn_attypes and at1 not in bkbn_attypes):
                    # find the corresponding row in atmat_sb
                    sb_row = atmat_sb.loc[(atmat_sb["attype1"]==at1) & (atmat_sb["attype2"]==at2)]
                    if len(sb_row) == 0:
                        sb_row = atmat_sb.loc[(atmat_sb["attype1"]==at2) & (atmat_sb["attype2"]==at1)]
                    if len(sb_row) > 0:
                        atmat.at[index, "exp_aver"] = sb_row["exp_aver"].values[0]
                        atmat.at[index, "probability"] = sb_row["probability"].values[0]

        # read sc-bkbn correction:
        if r2_corr_sb is not None:
            corr = pd.read_csv(os.path.join(_this_dir, r2_corr_sb), sep=r"\s+")
            for index, row in atmat.iterrows():
                atpair = row["atom_pair"]
                corr_row = corr.loc[(corr["atom_pair"]==atpair)]
                if len(corr_row) == 0:
                    corr_row = corr.loc[(corr["atom_pair"]==atpair)]
                if len(corr_row) > 0:
                    atmat.at[index, "probability"] = atmat.at[index, "probability"] - atmat.at[index, "probability"] *corr_row["r2"].values[0]
        return atmat

    def define_special_nonlocal_dict(self):
        """
        Define the special_nonlocal_dict based on the energy values in atmat.
        positive energy: attractive interaction, negative energy: repulsive interaction
        """
        special_nonlocal_dict = []
        for atom_pair in self.atmat["atom_pair"].unique():
            atp1, atp2 = atom_pair.split("_")
            # write only H-O interaction, skip all the rest
            if (atp1 == "H" and atp2 not in  H_interactions) or (atp2 == "H" and atp1 not in H_interactions):
                continue
            energy = self.atmat.loc[self.atmat["atom_pair"]==atom_pair, "energy"].values[0]
            sigma = self.atmat.loc[self.atmat["atom_pair"]==atom_pair, "exp_aver"].values[0] / 2**(1/6)
            if energy > 0:
                special_nonlocal_dict.append({
                    "atomtypes": ([atom_pair.split("_")[0]], [atom_pair.split("_")[1]]),
                    "interaction": "att",
                    "sigma": float(sigma),
                    "epsilon": float(energy),
                })
            else:
                special_nonlocal_dict.append({
                    "atomtypes": ([atom_pair.split("_")[0]], [atom_pair.split("_")[1]]),
                    "interaction": "rep",
                    "sigma": None,
                    "epsilon": float(self.c12_rep_df.loc[((self.c12_rep_df["atp1"]==atp1) & (self.c12_rep_df["atp2"]==atp2)) | ((self.c12_rep_df["atp1"]==atp2) & (self.c12_rep_df["atp2"]==atp1)), "c12"].values[0]),
                })
        return special_nonlocal_dict

    def nonlocal_matrix(self):
        """
        Calculate the energy values for the atomtype pairs based on the probability values and the p_threshold.
        """
        # sc-sc
        condition_sc_sc = ~(self.atmat["attype1"].isin(bkbn_attypes) | self.atmat["attype2"].isin(bkbn_attypes))
        self.atmat.loc[condition_sc_sc, "energy"] = np.array(self.regall(self.atmat.loc[condition_sc_sc, "probability"].to_numpy(), 1, self.emax, self.pth))
        # sc-bkbn
        condition_sc_bkbn = (self.atmat["attype1"].isin(bkbn_attypes) & ~self.atmat["attype2"].isin(bkbn_attypes)) | (self.atmat["attype2"].isin(bkbn_attypes) & ~self.atmat["attype1"].isin(bkbn_attypes))
        self.atmat.loc[condition_sc_bkbn, "energy"] = np.array(self.regall(self.atmat.loc[condition_sc_bkbn, "probability"].to_numpy(), 1, self.emax_sc_bkbn, self.pth_sc_bkbn))
        # bkbn-bkbn
        condition_bkbn_bkbn = self.atmat["attype1"].isin(bkbn_attypes) & self.atmat["attype2"].isin(bkbn_attypes)
        self.atmat.loc[ condition_bkbn_bkbn, "energy"] = np.array(self.regall(self.atmat.loc[condition_bkbn_bkbn, "probability"].to_numpy(), 1, self.emax_bkbn, self.pth_bkbn))

        self.special_nonlocal_dict = self.define_special_nonlocal_dict()
