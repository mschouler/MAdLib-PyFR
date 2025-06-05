import logging
import numpy as np
import os
import pandas as pd

from madlib_pyfr.adaptation_utils import get_files, run_cmd, time_function
from subprocess import run

logger = logging.getLogger(__name__)


class Adapter:
    """
    Mesh adaptation class.
    """
    def __init__(
            self,
            config: dict,
            step: int = 0,
            gmsh_file: str = "",
            sol_file: str = ""
    ):
        self.step = step
        self.config = config
        # the adaptation mode defines whether the adaptation
        # should be instantaneous (standard) or not (any other string)
        # Note: non-standard adaptation means that several solution timesteps
        #       are used by the adaptation script (i.e. metrics are combined)
        self.mode: str = config.get("mode", "standard")
        logger.info(f"adaptation mode set to: {self.mode}")
        self.outdir_pattern: str = config.get("outdir", "out")
        self.madlib_exe: str = config.get("madlibexe", "")
        self.vtk_exe: str = config.get("vtkexe", "")
        self.gmsh_file: str = gmsh_file if gmsh_file else config["msh_file"]
        self.restart: bool = True if sol_file else False

        # pyfr options
        pyfr = config["pyfr"].get("options", {})
        self.pyfr_opt: list[str] = [e for p in pyfr.items() for e in p]
        self.ini_file: str = config["pyfr"]["ini_file"]
        self.t: float = config["pyfr"].get("t0", 0.)
        self.tf: float = config["pyfr"]["tf"]
        # dt may be reduced if the simulation fails
        self.dt: float = config["pyfr"]["dt"][0]
        # pseudo_dt may be reduced if the simulation fails
        self.pseudo_dt: float = config["pyfr"].get("pseudo_dt", [-1])[0]
        # dt_sol is fixed
        self.dt_sol: float = config["pyfr"]["dt_sol"]
        # K is fixed unless the adaptive period control is used
        self.K: int = config["pyfr"].get("K", 1)
        # extra_dt
        self.extra_dt: float | None = config["pyfr"].get("extra_dt")
        self.extra_dt_sol: float | None = config["pyfr"].get("extra_dt_sol")
        # solution timestamp format
        self.format = config["pyfr"].get("format", ".2f")
        # convergence criteria
        self.ftt: int = config["pyfr"].get("ftt", 1)
        self.stat_crit: dict = config["pyfr"].get("statistical_criteria", {})
        self.QoI_list: list[np.ndarray] = []
        self.adaptation_list: list[int] = []

        # madlib options
        madlib = config["madlib"].get("options", {})
        self.madlib_opt: list[str] = [e for p in madlib.items() for e in p]
        self.timeout: float = config["madlib"].get("timeout", 5.)
        self.cmp_factor: float = config["madlib"].get("cmp_factor", 1)
        self.cmp: None | float = None

        self.sol_pattern: str = self.gmsh_file.split("/")[-1].split(".")[0]
        self.sol_file: str = sol_file if sol_file else self.get_sol_file()

        # files tracker
        self.dvtk: dict[int, str] = {}
        self.dmsh: dict[int, str] = {}
        self.dgmsh: dict[int, str] = {step: self.gmsh_file}
        self.dsol: dict[int, str] = {step: self.sol_file}
        self.derror: dict[int, float] = {}

        # timers
        self.t_run: float = 0.
        self.t_import: float = 0.
        self.t_export: float = 0.
        self.t_interpolate: float = 0.
        self.t_extvtk: float = 0.
        self.t_madlib: float = 0.

    def get_sol_file(self):
        """
        Returns the name of the solution file
        produced at the current adaptation step.
        """
        sol_file = os.path.join(
            self.get_outdir(),
            self.sol_pattern
            + f"_{self.step * self.K * self.dt_sol:{self.format}}.pyfrs"
        )
        return sol_file

    def get_outdir(self):
        """
        Returns the name of the output directory
        of the current adaptation step.
        """
        return self.outdir_pattern + f"_{self.step}"

    @time_function
    def import_mesh(self, gmsh_file: str) -> str:
        """
        Converts mesh from .msh to .pyfrm and
        returns the path to the new .pyfrm mesh.
        """
        outfile = gmsh_file.split("/")[-1][:-4] + ".pyfrm"
        pyfr_cmd = ["pyfr", "import",
                    gmsh_file, os.path.join(self.outdir, outfile)]
        run_cmd(pyfr_cmd, os.path.join(self.outdir, "import.out"))
        return os.path.join(self.outdir, outfile)

    def get_timestamp(self) -> str:
        """
        Returns the timestamp of the current adaptation step solution file.
        """
        return f"{self.t + self.K * self.dt_sol:{self.format}}"

    def update_ini(self, final_time: float) -> tuple[str, str]:
        """
        Creates new .ini file from the template and returns its path.
        """
        with open(self.ini_file, 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace(
            "@initial_time", str(self.t)
        )
        filedata = filedata.replace(
            "@final_time", str(final_time)
        )
        filedata = filedata.replace("@pseudo_dt", str(self.pseudo_dt))
        filedata = filedata.replace("@dt_sol", str(self.dt_sol))
        filedata = filedata.replace("@dt", str(self.dt))
        filedata = filedata.replace("@outdir", self.outdir)
        filedata = filedata.replace(
            "@sol_pattern", self.sol_pattern + f"_{{t:{self.format}}}"
        )
        # Write the file out again
        outfile = self.ini_file.split("/")[-1]
        with open(os.path.join(self.outdir, outfile), 'w') as file:
            file.write(filedata)
        return (
            os.path.join(self.outdir, outfile),
            os.path.join(self.outdir, self.sol_pattern
                         + f"_{self.get_timestamp()}.pyfrs")
        )

    @time_function
    def run_pyfr(self, mesh_file: str, ini_file: str, sol_file: str):
        """
        Runs pyfr from scratch or restarts from previous solution
        with the cuda backend and the given mesh and input files.
        For the first execution, the initial solution generated by pyfr
        is exported to .vtk.
        """
        if self.step == 1 and not self.restart:
            pyfr_cmd = ["pyfr", "run"] + self.pyfr_opt + [mesh_file, ini_file]
        else:
            pyfr_cmd = ["pyfr", "restart"] + self.pyfr_opt
            pyfr_cmd += [mesh_file, sol_file, ini_file]
        logger.debug(f"pyfr cmd: {pyfr_cmd}, step {self.step}")
        run_cmd(
            pyfr_cmd,
            os.path.join(self.outdir, f"run_{self.get_timestamp()}.out")
        )

    @time_function
    def export_vtk(self, mesh_file: str, sol_file: str = "") -> str:
        """
        Converts solutions from .pyfrs to .vtu and
        returns the path to the last .vtu file.
        """
        sol_list = [sol_file] if sol_file else get_files(self.outdir, ".pyfrs")
        for sol_file in sol_list:
            outfile = sol_file.split("/")[-1][:-6] + ".vtu"
            pyfr_cmd = ["pyfr", "export", "volume",
                        mesh_file, sol_file,
                        os.path.join(self.outdir, outfile)]
            run_cmd(pyfr_cmd, os.path.join(self.outdir, "export.out"))
        return os.path.join(self.outdir, outfile)

    @time_function
    def compute_vtk(self, vtk_file: str = "") -> str:
        """
        Computes VelocityMagnitude and VorticityMagnitude fields from
        the given .vtk solutions and
        returns the path to the last .vtu file.

        Note: in case of failure, already computed files are simply ignored.
        """
        vtk_list = (
            [vtk_file] if vtk_file else get_files(self.outdir, ".vtu", "-ext")
        )
        for file in vtk_list:
            ii = file.rfind("_")
            outfile = (file[:ii] + "-ext_" + file[ii + 1:]).split("/")[-1]
            vtk_cmd = [self.vtk_exe, file, os.path.join(self.outdir, outfile)]
            logger.debug(f"vtk_file: {file}, outfile: {outfile}")
            run_cmd(vtk_cmd, os.path.join(self.outdir, "vtk.out"))
            # this avoids duplicates as the method is called twice
            # on the same output directory
            if not vtk_file:
                self.dvtk[self.step] = os.path.join(self.outdir, outfile)
        return os.path.join(self.outdir, outfile)

    @time_function
    def adap_madlib(self, gmsh_file: str, txt_file: str) -> str:
        """
        Runs madlib with the given mesh (.msh) and solution (.vtu) files
        and returns the path to the adapted mesh (.msh).
        """
        t_stamp = self.get_timestamp()
        g_file = gmsh_file.split("/")[-1]
        ii = g_file.rfind("_") if "_" in g_file else g_file.rfind(".")
        outfile = g_file[:ii] + f"_{t_stamp}.msh"
        logger.debug(f"gmsh_file: {gmsh_file}, outfile: {outfile}")
        madlib_cmd = [self.madlib_exe, "-m", gmsh_file, "-f", txt_file,
                      "-o", os.path.join(self.outdir, outfile)]
        madlib_cmd.extend(self.madlib_opt)
        logger.debug(f"madlib cmd: {madlib_cmd}")
        stdfile = os.path.join(self.outdir, f"madlib_{t_stamp}.out")
        with open(stdfile, "wb") as out:
            run(madlib_cmd, stdout=out, stderr=out,
                check=True, timeout=self.timeout)
        return os.path.join(self.outdir, outfile)

    def edit_gmsh(self, gmsh_file: str, adapted_gmsh_file: str) -> str:
        """
        Changes the gmsh file format from 2.1 (produced with madlib)
        to 2.2 (required by pfr) and returns the path
        to the edited adapted gmsh file.

        Note: it also renumbers the nodes in a contiguous fashion.
        """
        gmsh_file_content = open(gmsh_file, "r").read().splitlines()
        adapted_gmsh_file_content = (
            open(adapted_gmsh_file, "r").read().splitlines()
        )
        # get original version
        v_index = gmsh_file_content.index("$MeshFormat") + 1
        # update version
        adapted_gmsh_file_content[v_index] = gmsh_file_content[v_index]
        # conserve physical names
        i_0 = gmsh_file_content.index("$PhysicalNames")
        i_n = gmsh_file_content.index("$EndPhysicalNames")
        adapted_gmsh_file_content = (
            adapted_gmsh_file_content[:i_0]
            + gmsh_file_content[i_0:i_n + 1]
            + adapted_gmsh_file_content[i_0:]
        )
        # make nodes number contiguous
        self.reorder_nodes(adapted_gmsh_file_content)
        # update adapted file content
        with open(adapted_gmsh_file, 'w') as ftw:
            ftw.write("\n".join(adapted_gmsh_file_content) + "\n")
        return adapted_gmsh_file

    def reorder_nodes(self, mesh_file_content: list[str]):
        """
        Makes the node index contiguous and update the elements acoordingly.

        Note: this is required by PyFR with point sampling.
        """
        # read mesh
        for idx, line in enumerate(mesh_file_content):
            if line.strip() == "$Nodes":
                node_idx = idx
            elif line.strip() == "$EndNodes":
                end_node_idx = idx
            elif line.strip() == "$Elements":
                elt_idx = idx
            elif line.strip() == "$EndElements":
                end_elt_idx = idx
        # get number of nodes and ids
        _ = int(mesh_file_content[node_idx + 1].strip())
        id_nodes = [int(n.strip().split()[0])
                    for n in mesh_file_content[node_idx + 2:end_node_idx]]
        # update node ids
        for idx, line in enumerate(mesh_file_content[node_idx + 2:end_node_idx]):
            new_line = " ".join([str(idx + 1)] + line.strip().split()[1:])
            mesh_file_content[node_idx + 2 + idx] = new_line
        # update element node ids
        for idx, line in enumerate(mesh_file_content[elt_idx + 2:end_elt_idx]):
            # nb node in the current line
            nb_node = (len(line.strip().split())
                       - (3 + int(line.strip().split()[2])))
            # compute new nodes ids
            new_nodes_id = [str(id_nodes.index(int(nid)) + 1)
                            for nid in line.strip().split()[-nb_node:]]
            # update line
            new_line = " ".join(line.strip().split()[:-nb_node] + new_nodes_id)
            mesh_file_content[elt_idx + 2 + idx] = new_line

    @time_function
    def resample_pyfr(
        self, mesh_file: str, sol_file: str,
        new_sol_file: str, adapted_mesh_file: str,
        ini_file: str
    ) -> str:
        """
        Interpolates the solution to the adapted mesh
        and returns the path to the interpolated solution.
        """
        outfile = new_sol_file.split("/")[-1]
        pyfr_cmd = [
            "pyfr", "resample", "-i", "idw",
            mesh_file, sol_file, adapted_mesh_file, ini_file, new_sol_file
        ]
        logger.debug(f"resample cmd: {pyfr_cmd}")
        with open(os.path.join(self.outdir, "resample.out"), "wb") as out:
            run(pyfr_cmd, stdout=out, stderr=out, check=True)
        return os.path.join(self.outdir, outfile)

    def get_txt_file(self, sol_file: str) -> str:
        """
        Writes and returns the file containing the list of solutions
        to be used in the adaptation.

        Note: in standard mode, only the last solution file is returned.
        """
        txt_file = os.path.join(self.outdir, f"list_{self.step}.txt")
        txt_content = (
            get_files(self.outdir, "-ext")[-self.K:]
            if not self.mode == "standard"
            else [sol_file]
        )
        with open(txt_file, "w") as tfile:
            tfile.write("\n".join(txt_content))
        return txt_file

    def get_error(self, outdir) -> float:
        """
        Reads and returns an estimate of the solution error.
        """
        txt = open(os.path.join(outdir, "error.txt"), "r").read().splitlines()
        return float(txt[0])

    def get_cmp(self, outdir) -> float:
        """
        Reads and returns the complexity of the mesh.
        """
        txt = open(os.path.join(outdir, "cmp.txt"), "r").read().splitlines()
        return float(txt[0])

    def adaptation_step(
            self,
            ext_vtk_file: str,
            new_sol_file: str,
            ini_file: str
    ):
        """
        Performs a complete adaptation iteration.

        Note: if ext_vtk_file, new_sol_file and inif_file are provided
              then the adaptation is performed right away.
        """
        if not (ext_vtk_file and new_sol_file and ini_file):
            self.step += 1

            self.outdir = self.get_outdir()
            os.makedirs(self.outdir, exist_ok=True)

            if self.step == 1:
                self.mesh_file, duration = self.import_mesh(self.gmsh_file)
                self.t_import += duration
                self.dmsh[self.step] = self.mesh_file
                logger.info(f"imported mesh_file: {self.mesh_file}")

            final_time = self.t + self.K * self.dt_sol
            ini_file, new_sol_file = self.update_ini(final_time)
            logger.info(f"ini_file: {ini_file}, new_sol_file: {new_sol_file}")

            # run pyfr
            _, duration = self.run_pyfr(
                self.mesh_file,
                ini_file,
                self.sol_file
            )
            self.t_run += duration
            logger.info(f"produced new_sol_file: {new_sol_file}")

            # export to vtk
            vtk_file, duration = self.export_vtk(self.mesh_file)
            self.t_export += duration
            logger.info(f"produced vtk_file: {vtk_file}")

            # compute derived quantities
            ext_vtk_file, duration = self.compute_vtk()
            self.t_extvtk += duration
            logger.info(f"produced ext_vtk_file: {ext_vtk_file}")

        # adaptation
        txt_file = self.get_txt_file(ext_vtk_file)
        adapted_gmsh_file, duration = self.adap_madlib(
            self.gmsh_file, txt_file
        )
        self.adaptation_list.append(self.step)
        self.t_madlib += duration
        logger.info(f"produced adapted_gmsh_file: {adapted_gmsh_file}")
        # read the complexity that will be enforced
        if self.cmp is None:
            self.cmp = self.get_cmp(self.outdir)
            logger.info(f"initial complexity: {self.cmp}")
            self.cmp *= self.cmp_factor
            logger.info(f"complexity set to: {self.cmp}")
            self.madlib_opt += ["--cmp", str(self.cmp)]
        else:
            self.cmp *= self.cmp_factor
            self.madlib_opt[-1] = str(self.cmp)
            logger.info(f"complexity set to: {self.cmp}")

        # edit and import mesh
        adapted_gmsh_file = self.edit_gmsh(self.gmsh_file, adapted_gmsh_file)
        logger.info(f"produced edited adapted_gmsh_file: {adapted_gmsh_file}")
        adapted_mesh_file, duration = self.import_mesh(adapted_gmsh_file)
        self.t_import += duration
        logger.info(f"produced adapted_mesh_file: {adapted_mesh_file}")

        # solution interpolation
        adapted_sol_file, duration = self.resample_pyfr(
            self.mesh_file,
            new_sol_file,
            new_sol_file,
            adapted_mesh_file,
            ini_file
        )
        self.t_interpolate += duration
        logger.info(f"produced adapted_sol_file: {adapted_sol_file}")

        # export to vtk
        # Note: the vtk solution is replaced by the solution interpolated
        #       on the adapted mesh
        vtk_file, duration = self.export_vtk(
            adapted_mesh_file, adapted_sol_file
        )
        self.t_export += duration
        logger.info(f"produced adapted vtk file: {vtk_file}")

        # compute derived quantities
        ext_vtk_file, duration = self.compute_vtk(vtk_file)
        self.t_extvtk += duration
        logger.info(f"produced ext_adapted vtk file: {ext_vtk_file}")

        # update inner attributes
        self.mesh_file = adapted_mesh_file
        self.dmsh[self.step] = adapted_mesh_file
        self.sol_file = adapted_sol_file
        self.dsol[self.step] = adapted_sol_file
        self.gmsh_file = adapted_gmsh_file
        self.dgmsh[self.step] = adapted_gmsh_file

        # increment time counter
        self.t += self.K * self.dt_sol
        logger.info(f"increment time counter t to: {self.t}")

    def decrement(self):
        """
        Decrements the adapter:
        - remove the gmsh mesh generated at the last step,
        - remove the pyfr solution generated at the last step,
        - remove the vtk solutions generated at the last step,
        - decrements the step counter.
        """
        _ = self.dgmsh.pop(self.step, None)
        logger.debug(f"removed from dgmsh: {_}")
        _ = self.dmsh.pop(self.step, None)
        logger.debug(f"removed from dmsh: {_}")
        _ = self.dsol.pop(self.step, None)
        logger.debug(f"removed from dsol: {_}")
        _ = self.dvtk.pop(self.step, None)
        logger.debug(f"removed from dvtk: {_}")
        self.step -= 1

    def extra_step(
            self, extra_dt: float, dt_sol: float
    ) -> tuple[str, str, str]:
        """
        Runs pyfr for an extra step without adaptation and returns
        the necessary files to make it possible to run an adaptation
        right after this step (useful for APC mode).
        """
        if extra_dt is None or dt_sol is None:
            logger.warning("extra_dt or extra_dt_sol not specified, "
                           "extra step aborted")
            return

        self.step += 1

        self.outdir = self.get_outdir()
        os.makedirs(self.outdir, exist_ok=True)

        final_time = self.t + extra_dt
        logger.info(f"final_time: {final_time}")
        self.dt_sol = dt_sol
        ini_file, new_sol_file = self.update_ini(final_time)
        logger.info(f"ini_file: {ini_file}, new_sol_file: {new_sol_file}")

        # run pyfr
        _, duration = self.run_pyfr(self.mesh_file, ini_file, self.sol_file)
        self.t_run += duration
        logger.info(f"produced new_sol_file: {new_sol_file}")

        # export to vtk
        vtk_file, duration = self.export_vtk(self.mesh_file)
        self.t_export += duration
        logger.info(f"produced vtk_file: {vtk_file}")

        # compute derived quantities
        ext_vtk_file, duration = self.compute_vtk()
        self.t_extvtk += duration
        logger.info(f"produced ext_vtk_file: {ext_vtk_file}")

        # update inner attributes
        self.dmsh[self.step] = self.mesh_file
        self.sol_file = new_sol_file
        self.dsol[self.step] = new_sol_file
        self.dgmsh[self.step] = self.gmsh_file

        # increment time counter
        self.t += extra_dt
        logger.info(f"increment time counter t to: {self.t}")

        return ext_vtk_file, new_sol_file, ini_file

    def statistical_phase(self) -> tuple[str, str, str]:
        """
        Runs pyfr to evacuate the transient phase resulting from
        the last adaptation step by monitoring the statistics of
        some quantities. Returns the necessary files to make it
        possible to run an adaptation right after this step
        (useful for APC mode).

        Note: in addition to the statistical criterion, we also make
              sure the phase is not run for more than one flow through time.
        """
        t0 = self.t
        max_t = self.ftt * self.stat_crit["max_ftt"]
        QoI_list: list[np.ndarray] = []
        while True:
            # compute extra step
            ext_vtk_file, new_sol_file, ini_file = self.extra_step(
                self.K * self.dt_sol,
                self.dt_sol
            )
            # compute monitored values
            QoI_sublist = []
            for f, k in zip(self.stat_crit["files"], self.stat_crit["keys"]):
                df = pd.read_csv(f)
                # only compute the mean over the solutions in the last step
                QoI_sublist.append(np.nanmean(
                    df.loc[df["t"] >= t0][k],
                    axis=0
                ))
            QoI_list.append(np.concatenate(QoI_sublist))
            # compute error
            err = self.stat_cv(QoI_list)
            self.derror[self.step] = err
            logger.info(f"transient convergence criteria: {err}")
            # check termination criteria
            if err < self.stat_crit["criteria"]:
                logger.info("statistical phase converged")
                break
            if self.t > t0 + max_t * self.dt:
                logger.info("statistical phase reached its time limit")
                break

        # save the converged QoIs of this adaptation step
        logger.info(f"statistical phase finished with QoI_list: {QoI_list}")
        self.QoI_list.append(QoI_list[-1])

        return ext_vtk_file, new_sol_file, ini_file

    def convergence(self) -> bool:
        """
        Returns True if the mesh adaptation should stop, False otherwise.
        """
        # basic case
        if not self.mode == "APC":
            return self.t >= self.tf
        # automatic period control
        else:
            err = self.stat_cv(self.QoI_list)
            logger.info(f"global convergence criteria: {err}")
            return err < self.stat_crit["criteria"]

    def stat_cv(self, y) -> float:
        """
        Returns the error criterion.

        Note: the convergence criteria is a percentage computed
              over three iterations.
        """
        if len(y) > 2:
            logger.debug(f"y in stat_cv: {y}")
            logger.debug(f"y1 - y2 / y1: {abs((y[-1] - y[-2]) / y[-1])}")
            logger.debug(f"y1 - y3 / y1: {abs((y[-1] - y[-3]) / y[-1])}")
            err = max(
                np.max(abs((y[-1] - y[-2]) / y[-1])),
                np.max(abs((y[-1] - y[-3]) / y[-1]))
            ) * 100
            return err
        else:
            return np.nan

    def adjust_controller(self):
        """
        Reduce or extend K depending on the error.
        """
        qoi = self.QoI_list
        logger.info(f"qoi: {qoi}")
        if len(qoi) < 2:
            return
        epsilon = np.min((qoi[-1] - qoi[-2]) / (qoi[-1] + qoi[-2]))
        logger.info(f"APC criterion: {epsilon}")
        # self.K = int((1 - epsilon) * self.K)
        logger.info(f"APC adapted K to new value: {self.K}")
