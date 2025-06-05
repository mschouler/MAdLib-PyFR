import argparse
import logging
import os
import time

from madlib_pyfr.adaptation_utils import get_files, set_logger, time_function, write_series
from madlib_pyfr.adaptation import Adapter
from subprocess import run

logger = logging.getLogger(__name__)

MADLIB = "/home/mschouler/Documents/Sorbonne/MAdLib_2_0_0/trunk/build/Testcases/2dEulerVortex"  # noqa


class EulerAdapter(Adapter):
    def __init__(self, mesh_file: str, ini_file: str, nite: int, dt_sol: float, outdir: str):
        # build dummy config
        config = {
            "mode": "standard",
            "msh_file": mesh_file,
            "outdir": outdir,
            "pyfr": {
                "options": {
                    "-b": "openmp"
                },
                "ini_file": ini_file,
                "t0": 0.,
                "dt": [dt_sol],
                "dt_sol": dt_sol,
                "tf": dt_sol * nite,
                "format": ".0f"
            },
            "madlibexe": MADLIB,
            "madlib": {}
        }
        super().__init__(config)

    def get_cmp(self, outdir) -> float:
        """
        Reads and returns the complexity of the mesh.

        Note: for this is useless.
        """
        return -1

    def get_txt_file(self, sol_file: str) -> str:
        """
        Writes and returns the file containing the list of solutions
        to be used in the adaptation.

        Note: for this example the solution file is directly passed to MAdLib.
        """
        return sol_file

    @time_function
    def adap_madlib(self, gmsh_file: str, sol_file: str) -> str:
        """
        Runs madlib with the given mesh (.msh) and solution (.vtu) files
        and returns the path to the adapted mesh (.msh).

        Note: for this example the adaptation script only takes two positional input arguments.
        """
        t_stamp = self.get_timestamp()
        g_file = gmsh_file.split("/")[-1]
        ii = g_file.rfind("_") if "_" in g_file else g_file.rfind(".")
        outfile = g_file[:ii] + f"_{t_stamp}.msh"
        logger.debug(f"gmsh_file: {gmsh_file}, outfile: {outfile}")
        madlib_cmd = [self.madlib_exe, gmsh_file, sol_file, os.path.join(self.outdir, outfile)]
        logger.debug(f"madlib cmd: {madlib_cmd}")
        stdfile = os.path.join(self.outdir, f"madlib_{t_stamp}.out")
        with open(stdfile, "wb") as out:
            run(madlib_cmd, stdout=out, stderr=out,
                check=True, timeout=self.timeout)
        return os.path.join(self.outdir, outfile)

    @time_function
    def compute_vtk(self, vtk_file: str = "") -> str:
        """
        Computes VelocityMagnitude and VorticityMagnitude fields from
        the given .vtk solutions and
        returns the path to the last .vtu file.

        Note: in this example, no such step is necessary.
        """
        vtk_file_list = get_files(self.outdir, ".vtu", "-ext")
        logger.debug(f"vtk_file_list: {vtk_file_list}")
        # the first step contains two solutions: 0.vtu and 1.vtu
        for i, file in enumerate(vtk_file_list[::-1]):
            self.dvtk[self.step - i] = file
        return vtk_file_list[-1]


def main():
    """
    Couples PyFR with MAdLib.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--mesh", type=str, help="gmsh file", required=True)
    parser.add_argument("-i", "--ini", type=str, help="PyFR input file", required=True)
    parser.add_argument("-n", "--nite", type=int, help="the number of iterations", default=10)
    parser.add_argument("-out", "--outdir", type=str, help="output directory", default="output/out")
    parser.add_argument(
        "-dt", "--time-step", type=float, help="the saved solution timestep", default=1.)
    args = parser.parse_args()

    # set logger
    log_level = getattr(logging, "DEBUG")
    set_logger("adap.log", log_level)

    t0 = time.time()

    # initialize Adapter
    mesh_adapter = EulerAdapter(args.mesh, args.ini, args.nite, args.time_step, args.outdir)

    # adaptation loop
    vtk_file, sol_file, ini_file = "", "", ""
    for ii in range(args.nite):
        logger.info(f"\nITERATION {mesh_adapter.step + 1}..")
        mesh_adapter.adaptation_step(vtk_file, sol_file, ini_file)

    adap_duration = time.time() - t0
    logger.info(f"\nadaptation finished in {adap_duration:.2f} sec.")

    # create .series file for paraview
    series_file = write_series(mesh_adapter.sol_pattern, mesh_adapter.dvtk, os.getcwd())
    logger.info(f"\nproduced series file: {series_file}")

    # write execution time statistics
    logger.info("\n\n **EXECUTION TIME STATISTICS**")
    logger.info(f"adaptation time: {adap_duration:.2f} sec. (100 %)")
    logger.info(f"pyfr execution time: {mesh_adapter.t_run:.2f} sec. "
                f"({mesh_adapter.t_run / adap_duration * 100:.2f} %)")
    logger.info(f"pyfr import time: {mesh_adapter.t_import} sec. "
                f"({mesh_adapter.t_import / adap_duration * 100:.2f} %)")
    logger.info(f"pyfr export time: {mesh_adapter.t_export:.2f} sec. "
                f"({mesh_adapter.t_export / adap_duration * 100:.2f} %)")
    logger.info(f"pyfr interpolation time: {mesh_adapter.t_interpolate:.2f}"
                f"sec. "
                f"({mesh_adapter.t_interpolate / adap_duration * 100:.2f} %)")
    logger.info(f"vtk extension time: {mesh_adapter.t_extvtk:.2f} sec. "
                f"({mesh_adapter.t_extvtk / adap_duration * 100:.2f} %)")
    logger.info(f"madlib time: {mesh_adapter.t_madlib:.2f} sec. "
                f"({mesh_adapter.t_madlib / adap_duration * 100:.2f} %)")


if __name__ == '__main__':
    main()
