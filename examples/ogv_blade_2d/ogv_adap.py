import logging
import numpy as np
import os
import pandas as pd
import sys

from madlib_pyfr.adaptation import Adapter

logger = logging.getLogger(__name__)

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from inc_cylinder_2d.cylinder_adap import main as base_main


def compute_pts(nb_pts: int, boundaries: list[list[float]]) -> str:
    """
    Returns a list of nb_pts evenly spaced points along the line between the boundaries.

    Note: the returned value is a string '[(x0, y0, z0), ..., (x_nb_pts, y_nb_pts, z_nb_pts)]'
          to be fed to PyFR plugin sampler 
    """
    x = np.array(boundaries[0])
    y = np.array(boundaries[1])  
    t = np.linspace(0, 1, nb_pts + 1, endpoint=False)[1:]
    XY = x + t[:, np.newaxis] * (y - x)
    return str(list(map(tuple, XY.tolist())))


def compute_mixedout_qty(df: pd.DataFrame) -> tuple[float, float]:
    """
    Returns the mixed-out static and total pressures:
    see A. Prasad (2004): https://doi.org/10.1115/1.1928289
    """
    # conservation of mass
    m_bar = np.nanmean(df["u"] * df["rho"])
    rho_bar = np.nanmean(df["rho"])
    p_bar = np.nanmean(df["p"])
    u_bar = m_bar / rho_bar
    v_bar = np.nanmean(df["v"]) / rho_bar
    uu_bar = u_bar**2
    vv_bar = v_bar**2

    # conservation of momentum
    x_mom = m_bar * u_bar + p_bar
    y_mom = m_bar * v_bar

    # conservation of energy
    gamma = 1.4
    R = 287.058
    E = m_bar * gamma / (gamma - 1) * p_bar / rho_bar + m_bar / 2. * (uu_bar + vv_bar)

    # quadratic equation
    Q = 1 / m_bar**2 * (1 - 2 * gamma / (gamma - 1))
    L = 2 / m_bar**2 * (gamma / (gamma - 1) * x_mom - x_mom)
    C = 1 / m_bar**2 * (x_mom**2 + y_mom**2) - 2 * E / m_bar

    # select subsonic root
    p_bar = (-L - np.sqrt(L**2 - 4 * Q * C)) / 2 / Q
    T_bar = p_bar / rho_bar / R
    T0_bar = (gamma - 1) / (gamma * R) * E / m_bar
    p0_bar = p_bar * (T0_bar / T_bar)**(gamma / (gamma - 1))
    return p_bar, p0_bar


class OGVAdapter(Adapter):
    def __init__(self, *args, **kwargs):
        """
        This class enhances the Adapter class to allow for the computation
        of the mixed-out loss coefficient of the cascade.
        """
        print("OGVAdapter init")
        super().__init__(*args, **kwargs)
        # build interpolation lines
        self.interpolation_lines: list[list[float]] = self.config["pyfr"]["interpolation_lines"]
        nb_pts = self.interpolation_lines["nb_pts"]
        input_line = self.interpolation_lines["input_line"]
        output_line = self.interpolation_lines["output_line"]
        self.input_pts = compute_pts(nb_pts, input_line)
        self.output_pts = compute_pts(nb_pts, output_line)
    
    def update_ini(self, final_time: float, filedata: str = "") -> tuple[str, str]:
        """
        Creates new .ini file from the template and returns its path.

        Note: for this use-case, the sampling points are computed based on the user specifications.
        """
        with open(self.ini_file, 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace(
            "@input", self.input_pts
        )
        filedata = filedata.replace(
            "@output", self.output_pts
        )
        return super().update_ini(final_time, filedata)

    def compute_QoI_sublist(self, t0: float) -> np.ndarray:
        """
        Returns the mixed-out loss coefficient of the current step.
        """
        QoI_sublist = []
        input_df = pd.read_csv("input.csv")
        output_df = pd.read_csv("output.csv")
        # compute mixed-out quantities
        p1, p01 = compute_mixedout_qty(input_df.loc[input_df["t"] >= t0])
        _, p02 = compute_mixedout_qty(output_df.loc[output_df["t"] >= t0])
        # compute mixed-out loss
        QoI_sublist.append((p01 - p02) / (p01 - p1))
        return np.atleast_1d(QoI_sublist)

    def adaptation_step(self, ext_vtk_file, new_sol_file, ini_file):
        # original adaptation step
        super().adaptation_step(ext_vtk_file, new_sol_file, ini_file)
        # compute mixed-out loss coefficient if standard computation
        w = self.compute_QoI_sublist(self.t0)[-1]
        logger.info(f"step {self.step}, averaged mixed-out loss coefficient: {w}")
        # update QoI_list if necessary
        if self.mode != "APC":
            self.QoI_list.append(w)


def main():
    """
    Couples PyFR with MAdLib.

    Note: this main is just a wrapper around the cylinder coupling main.
    """
    base_main(Adapter)


if __name__ == '__main__':
    main()
