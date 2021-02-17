import gromologist as gml
from copy import deepcopy
import os
from subprocess import call
from multiprocessing import Pool

try:
    from scipy.optimize import dual_annealing
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    pass


class DihOpt:
    def __init__(self, top, qm_ref, cutoff=None, qm_unit='kj/mol', processes=1, tmpi=False, **kwargs):
        """
        The DihOpt is a dihedral optimizer that reads in (preferably) a Gaussian .log
        file and a Gromacs topology with selected dihedral terms marked as DIHOPT in the comment field,
        and then optimizes the dihedral parameters to minimize the RMSE between the QM and MM energies.
        :param top: a gml.Top object or str, Gromacs topology (either read-in or filename)
        :param qm_ref: str or list of str, path(s) to Gaussian .log files
        :param cutoff: float, will not attempt to fit data points above this threshold
        :param qm_unit: str, if qm_ref is a list of energies, this should be 'kcal/mol', 'kj/mol' or 'hartree'
        :param processes: int, how many independent fitting procedures should be launched (supports parallelism)
        :param tmpi: bool, whether to use Gromacs with thread-MPI in parallel execution
        :param kwargs: dict, if 'top' is a path, this will be passed to the gml.Top constructor
        """
        self.orig_top = top if isinstance(top, gml.Top) else gml.Top(top, **kwargs)
        self.orig_vals = None
        self.best_rmse = None
        self.best_top = self.orig_top
        self.optimizable_params = self.orig_top.parameters.get_opt_dih()
        self.traj = None
        self.frame = None
        self.lowest_iter_rmse = 1000000
        self.curr_iter = 0
        self.opt_results = None
        self.print_progress = False
        self.all_opts = None
        self.energy_term = 0
        self.energy_profiles_opt = []
        self.check_imports()
        # TODO so far only works with Gaussian outputs
        qm_unit = qm_unit.lower()
        if qm_unit not in ['kj/mol', 'kcal/mol', 'hartree']:
            raise RuntimeError("Unit can be kj/mol, kcal/mol or hartree")
        # we're converting everything to kJ/mol
        scaling = 1.0 if qm_unit == 'kj/mol' else 4.184 if 'kcal/mol' else 2625.5
        self.cutoff = cutoff * scaling if cutoff is not None else 10000000
        if isinstance(qm_ref, str) and qm_ref.endswith('.dat'):  # a file with energy values
            self.qm_ref = [float(x.split()[0])*scaling for x in open(qm_ref)]
        elif isinstance(qm_ref, str) and qm_ref.endswith('.log'):  # a .log file from Gaussian
            self.qm_ref = self.read_gaussian_energies(qm_ref)
        elif isinstance(qm_ref, list) and type(qm_ref[0]) in [int, float]:  # a list with energy values
            self.qm_ref = [x*scaling for x in qm_ref]
        elif isinstance(qm_ref, list) and type(qm_ref[0]) == str:  # a list of Gaussian logs
            self.qm_ref = [e for eners in [self.read_gaussian_energies(ref_file) for ref_file in qm_ref] for e in eners]
        else:
            raise RuntimeError("For qm_ref, pass a Gaussian .log file, a list of Gaussian .logs, a .dat file with "
                               "per-frame energies, or a list with energy values.")
        self.qm_ref = [x - min(self.qm_ref) for x in self.qm_ref]
        self.mod_tops = [deepcopy(self.orig_top) for _ in range(processes)]
        self.processes = processes
        self.run_count = 0
        self.tmpi = tmpi
        self.gmx = self.find_gmx()
        self.gen_dirs()
        self.write_mdp()
        self.opt_prof = None

    def check_imports(self):
        """
        To keep the lib as lightweight as possible, we're not explicitly requiring scipy/numpy,
        but to use DihOpt both are required, so here we look for them before running
        :return: None
        """
        try:
            _ = dual_annealing
        except NameError:
            raise ImportError("scipy needed to run dihedral optimization, run 'pip install scipy' in command line")
        try:
            _ = np.array
        except NameError:
            raise ImportError("numpy needed to run dihedral optimization, run 'pip install numpy' in command line")

    def calc_energy(self, ff_values=None, sys=0, cleanup=True):
        """
        Calculates the energy profile for a given set of dihedral parameters
        by running a rerun and extracting data from the resulting .xvg file
        :param ff_values: list of floats, dihedral parameters that will be set for the current iteration
        :param sys: int, index of the system
        :param cleanup: bool, whether to remove all files after the iteration
        :return: np.array, classical (MM) potential energies in kJ/mol for each frame in the trajectory
        """
        ff_values = self.best_top.parameters.get_opt_dih() if ff_values is None else ff_values
        self.mod_tops[sys].parameters.set_opt_dih(ff_values)
        self.mod_tops[sys].save_top('opt{}/mod.top'.format(sys))
        self.gmx_grompp(sys)
        self.gmx_mdrun(sys)
        energy = self.gmx_energy(sys)
        if cleanup:
            self.cleanup(sys)
        return energy

    def calc_rmse(self, ff_values=None, sys=0, cleanup=True):
        """
        For a given set of dihedral parameters, runs the rerun and calculates the RMSE
        between the QM data and the new MM energies
        :param ff_values: list of floats, dihedral parameters that will be set for the current iteration
        :param sys: int, index of the system
        :param cleanup: bool, whether to remove all files after the iteration
        :return: float, root-mean-square-error (in kJ/mol) between the QM and MM potential energy profiles
        """
        ff_values = self.best_top.parameters.get_opt_dih() if ff_values is None else ff_values
        energies = self.calc_energy(ff_values, sys, cleanup)
        rmse = [(x-y)**2 for x, y in zip(energies, self.qm_ref) if y < self.cutoff]
        rmse = (sum(rmse)/len(rmse))**0.5
        call('echo {r} >> opt{s}/rmse'.format(s=sys, r=rmse), shell=True)
        if self.curr_iter % 10 == 0 and int(self.curr_iter) > 0 and self.print_progress:
            print("Iteration {} in process {}, current lowest RMSE is {}".format(int(self.curr_iter), sys,
                                                                                 self.lowest_iter_rmse))
            self.print_progress = False
        if abs(self.curr_iter % 10) > 0 and not self.print_progress:
            self.print_progress = True
        return rmse

    def get_bounds(self, sys=0):
        """
        Finds bounds for parameter optimization, odd ones are for angle offsets,
        even ones for force constants
        :param sys: int, index of the system
        :return: list of tuples, bounds for the optimization procedure
        """
        vals = self.mod_tops[sys].parameters.get_opt_dih()
        return [(0, 360), (0, 15)] * (len(vals) // 2)

    def optimize(self, maxiter=200):
        """
        Calculates initial values, then performs the actual optimization
        using sklearn.optimize.dual_annealing either in a single process
        (when self.processes is 1) or in parallel (when self.processes > 1)
        :param maxiter: int, at how many iterations to terminate
        :return: None
        """
        # check for imports first
        self.run_count += 1
        self.orig_vals = self.calc_energy() if self.orig_vals is None else self.orig_vals
        self.best_rmse = self.calc_rmse()
        print("Running optimization in {} process{}, initial RMSE "
              "is {}".format(self.processes, 'es' if self.processes > 1 else '', self.best_rmse))
        if self.processes == 1:
            self.opt_results = dual_annealing(self.calc_rmse, bounds=self.get_bounds(0), maxiter=maxiter,
                                              x0=self.best_top.parameters.get_opt_dih(), callback=self.progress)
            # opt_vals = shgo(self.calc_rmse, args=(0,), bounds=self.get_bounds(0), options={'maxev': 1000})
            self.best_top = self.mod_tops[0]
            print("\nOptimized parameters:\n")
            print("    RMSE:       {:10.3f}".format(self.opt_results.fun))
            print("    Parameters: {}\n".format([round(x, 6) for x in self.opt_results.x]))
            print("Initial parameters:\n")
            print("    RMSE:       {:10.3f}".format(self.best_rmse))
            print("    Parameters: {}\n".format([round(x, 6) for x in self.orig_top.parameters.get_opt_dih()]))
        else:
            p = Pool()
            self.all_opts = p.map(mappable, [(i, self, maxiter) for i in range(self.processes)])
            best_model_index = int(np.argmin([q.fun for q in self.all_opts]))
            self.opt_results = self.all_opts[best_model_index]
            self.best_top = self.mod_tops[best_model_index]
            print("\nObtained {} models with the following RMSE values:\n".format(len(self.all_opts)))
            for i in range(len(self.all_opts)):
                print("    Model {}".format(i))
                print("    RMSE:       {:10.3f}".format(self.all_opts[i].fun))
                print("    Parameters: {}\n".format([round(x, 6) for x in self.all_opts[i].x]))
            print("Initial parameters:\n")
            print("    RMSE:       {:10.3f}".format(self.best_rmse))
            print("    Parameters: {}\n".format([round(x, 6) for x in self.orig_top.parameters.get_opt_dih()]))
            print("Choosing the lowest-RMSE model {}".format(best_model_index))
        self.best_top.parameters.set_opt_dih(self.opt_results.x)
        self.best_top.save_top('opt{}_topol.top'.format(self.run_count))
        self.energy_profiles_opt.append(self.calc_energy(self.opt_results.x))

    def restart(self, maxiter=200):
        """
        Can be used to reoptimize the parameters, using previous
        best result as the new starting point
        :param maxiter: int, at how many iterations to terminate
        :return:
        """
        self.best_top = gml.Top('opt{}_topol.top'.format(self.run_count))
        self.optimize(maxiter=maxiter)

    def gmx_grompp(self, sys=0):
        """
        Runs gmx grompp to prepare a .tpr file for the rerun
        :param sys: int, index of the system
        :return:
        """
        call('{gmx} grompp -f run.mdp -p opt{sys}/mod.top -c {frame} -o opt{sys}/mod.tpr -po opt{sys}/mdout.mdp '
             '-maxwarn 5 >> gmp.log 2>&1'.format(gmx=self.gmx, sys=sys, frame=self.frame), shell=True)

    def gmx_mdrun(self, sys=0):
        """
        Runs gmx mdrun with the '-rerun' option to extract MM energies
        :param sys: int, index of the system
        :return:
        """
        if self.tmpi:
            tmpi = ' -ntmpi 1 '
        else:
            tmpi = ''
        call('{gmx} mdrun -s opt{sys}/mod.tpr -o opt{sys}/mod.trr -x opt{sys}/mod.xtc -e opt{sys}/mod.edr '
             '-g opt{sys}/mod.log -c opt{sys}/mod.gro -v -ntomp 1 {tmpi} -rerun {traj} >> mdr.log '
             '2>&1'.format(gmx=self.gmx, sys=sys, tmpi=tmpi, traj=self.traj), shell=True)

    def gmx_energy(self, sys=0):
        """
        Runs gmx energy on the .edr file produced by the rerun, also
        extracts energies (in kJ/mol) from it & shifts them to 0
        :param sys: int, index of the system
        :return: numpy.array, energies of individual frames in kJ/mol
        """
        if not self.energy_term:
            self.get_energy_term(sys)
        call('echo "{et} 0" | {gmx} energy -f opt{sys}/mod.edr -o opt{sys}/mod.xvg >> edr.log '
             '2>&1'.format(gmx=self.gmx, sys=sys, et=self.energy_term), shell=True)
        vals = np.loadtxt('opt{sys}/mod.xvg'.format(sys=sys), comments=['#', '@'])[:, 1]
        return vals - np.min(vals)

    def get_energy_term(self, sys=0):
        """
        Checks the index of the 'Potential' term in gmx energy
        (usually 8 or 9, depending on the presence of impropers)
        :param sys: int, index of the system
        :return: None
        """
        call('echo "0" | {gmx} energy -f opt{sys}/mod.edr > legend.log 2>&1'.format(gmx=self.gmx, sys=sys), shell=True)
        self.energy_term = os.popen('cat legend.log | grep Potential | awk \'{for (i=1;i<=NF;i++){if ($i=="Potential")'
                                    ' {print $(i-1)}}}\'').read().strip()

    @staticmethod
    def cleanup(sys=0):
        """
        Removes files on-the-fly to prevent errors at 99 backups
        :param sys: int, index of the system
        :return: None
        """
        call('rm opt{s}/mdout.mdp opt{s}/mod.xvg opt{s}/mod.edr opt{s}/mod.tpr opt{s}/mod.log'.format(s=sys),
             shell=True)

    def plot_fit(self):
        """
        If optimization is complete, plots the resulting energy profile against
        both original MM energies and the QM reference
        :return: None
        """
        import matplotlib.pyplot as plt
        if self.opt_results is not None:
            plt.clf()
            plt.plot(self.orig_vals, label='original MM')
            for n, array in enumerate(self.energy_profiles_opt, 1):
                plt.plot(array, label='optimized MM (round {})'.format(n))
            plt.plot(self.qm_ref, label='QM reference', c='k')
            plt.legend()
            plt.show()

    def progress(self, _, f, context):
        """
        Constantly updates the curr_iter and lowest_iter_rmse attributes used for reporting
        :param _: discarded
        :param f: float, current function value (RMSE)
        :param context: int, key corresponding to different outcomes in the optimization algorithm
        :return: None
        """
        if context in [0, 1, 2]:
            self.curr_iter += 1.0
            if f < self.lowest_iter_rmse:
                self.lowest_iter_rmse = f

    @staticmethod
    def write_mdp():
        """
        Writes an .mdp file with all options necessary for a rerun
        :return: None
        """
        mdp_defaults = {"nstenergy": 1, "nstlog": 0, "nstcalcenergy": 1, "nstxout-compressed": 0, "pbc": "xyz",
                        "coulombtype": "Cut-off", "vdw-type": "Cut-off", "rlist": 2.0, "rcoulomb": 2.0, "rvdw": 2.0}
        mdp = '\n'.join(["{} = {}".format(param, value) for param, value in mdp_defaults.items()])
        with open('run.mdp', 'w') as outfile:
            outfile.write(mdp)

    def gen_dirs(self):
        """
        Makes separate directories for all the subprocesses
        :return: None
        """
        for i in range(self.processes):
            try:
                os.mkdir('opt{}'.format(i))
            except FileExistsError:
                pass

    @staticmethod
    def find_gmx():
        """
        Attempts to find Gromacs executables to run gmx grompp/mdrun/energy
        :return: str, path to the gmx/gmx_mpi/gmx_d executable
        """
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
        return gmx

    def read_gaussian_energies(self, log, traj='opt_traj.pdb', structfile='struct.pdb'):
        """
        Function that parses a Gaussian .log file with a dihedral scan results,
        extracting optimized geometries and reference QM energies; also saves
        the trajectory/structure as .pdb files based on .top atom names
        :param log: str, name of the Gaussian .log file
        :param traj: str, name of the output .pdb trajectory file
        :param structfile: str, name of the output .pdb structure file needed for gmx grompp
        :return: list, energies in kJ/mol
        """
        log_contents = [x.strip().split() for x in open(log)]
        natoms = [int(x[1]) for x in log_contents if x and x[0] == "NAtoms="][0]
        energies = []
        read_geo = False
        geoms = []
        energy_tmp = None
        for n, line in enumerate(log_contents):
            if len(line) > 2 and line[0] == 'SCF' and line[1] == 'Done:':
                energy_tmp = float(line[4])
            if len(line) > 3 and line[1] == 'Stationary' and line[2] == 'point' and line[3] == 'found.':
                if energy_tmp is not None:
                    energies.append(energy_tmp)
                else:
                    raise RuntimeError("Keyword 'Stationary point found' before first 'SCF done' in Gaussian log, check"
                                       " your QM run")
                read_geo = True
            if read_geo and len(line) > 1 and line[0] == 'Standard' and line[1] == 'orientation:':
                struct = np.array([[float(line[x]) for x in [3, 4, 5]] for line in log_contents[n+5:n+5+natoms]])
                geoms.append(struct)
                read_geo = False
        opt_traj = gml.Traj(top=self.orig_top, array=geoms)
        opt_traj.save_traj_as_pdb(traj)
        opt_traj.save_pdb(structfile)
        self.traj = traj
        self.frame = structfile
        return [2625.5 * x for x in energies]

    def make_movie(self):
        """
        Uses Molywood to produce a movie showing the optimization results
        along with the QM optimization trajectory, highlighting the dihedral
        being optimized
        :return: None
        """
        dihedrals_indices = self.orig_top.parameters.get_opt_dih_indices()
        np.savetxt('energies.dat', np.vstack([np.linspace(0, len(self.orig_vals), len(self.orig_vals))] + [self.qm_ref]
                                             + [self.orig_vals] + self.energy_profiles_opt).T)
        moly_inp = "$ global fps={} name=dihopt\n$ scene structure={} resolution=700,700" \
                   "\n\n#scene\n".format(len(self.orig_vals)/8, self.traj)
        moly_inp = moly_inp + 'highlight selection="all" style=licorice color=type mode=u\n'
        for dih in dihedrals_indices:
            moly_inp = moly_inp + 'highlight selection="serial {} {} {} {}" style=licorice thickness=1.1 ' \
                                  'alpha=0.65 mode=u color=green\n'.format(*dih)
        moly_inp = moly_inp + 'insert_tcl code="display resetview"\nzoom_out scale=1.5\nanimate frames=0\n'
        moly_inp = moly_inp + 'do_nothing t=2s\n'
        moly_inp = moly_inp + '{{animate t=8s frames=0:{}; add_overlay datafile=energies.dat origin=0.65,0.65 ' \
                              'relative_size=0.35}}\n'.format(len(self.orig_vals))
        moly_inp = moly_inp + 'do_nothing t=2s\n'
        with open('moly.inp', 'w') as outfile:
            outfile.write(moly_inp)
        call('molywood moly.inp', shell=True)


def mappable(arg):
    """
    Wrapper for dual_annealing needed to make it parallelizable
    through Pool.map() in DihOpt.optimize() (needs to be defined in the
    global scope to be picklable)
    :param arg: tuple, combines system index, the DihOpt object, and the maxiter value
    :return: result of the dual_annealing fn
    """
    sys, dihopt, maxiter = arg
    return dual_annealing(dihopt.calc_rmse, args=(sys, True), bounds=dihopt.get_bounds(0), maxiter=maxiter, seed=sys,
                          x0=dihopt.best_top.parameters.get_opt_dih(), callback=dihopt.progress)
