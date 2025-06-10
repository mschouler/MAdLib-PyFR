import argparse
import json
import logging
import os
import pandas as pd
import time
import traceback

from madlib_pyfr.adaptation_utils import set_logger, write_series
from madlib_pyfr.adaptation import Adapter

logger = logging.getLogger(__name__)


def main():
    """
    Couples PyFR with MAdLib.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--config", type=str, required=True, help="config file")
    parser.add_argument("--mesh", type=str, default="", help="gmsh file in case of restart")
    parser.add_argument(
        "--restart", type=int, default=0, help="restart adaptation from the given adaptation step")
    parser.add_argument(
        "--sol", type=str, default="", help=".pyfrs solution file in case of restart")
    parser.add_argument(
        "--verbose", type=str, default="DEBUG",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level")
    args = parser.parse_args()

    # set logger
    log_level = getattr(logging, args.verbose.upper())
    set_logger("adap.log", log_level)

    t0 = time.time()

    # load and edit config
    with open(args.config) as jfile:
        config = json.load(jfile)

    # initialize Adapter
    mesh_adapter = Adapter(config, step=args.restart, gmsh_file=args.mesh, sol_file=args.sol)

    # adaptation loop
    nfail = 0
    abort = False
    vtk_file, sol_file, ini_file = "", "", ""
    while not mesh_adapter.convergence():
        logger.info(f"\nITERATION {mesh_adapter.step + 1}..")
        try:
            # adaptation step
            mesh_adapter.adaptation_step(vtk_file, sol_file, ini_file)
            if mesh_adapter.mode == "APC":
                # statistical phase
                vtk_file, sol_file, ini_file = mesh_adapter.statistical_phase()

                # decrement time counter
                # Note: this ensures the solutions produced in the last
                #       step are re-used for the next adaptation step
                mesh_adapter.t -= mesh_adapter.K * mesh_adapter.dt_sol
                logger.info(f"decrement time counter t to: {mesh_adapter.t}")

                # adjust adaptation controller
                mesh_adapter.adjust_controller()

        except Exception as e:
            logger.error(f"adaptation failed with exception: {e}")
            logger.error(f"Traceback message:\n{traceback.format_exc()}")
            nfail += 1
            mesh_adapter.decrement()

            # fault tolerance protocol
            # Note: this was designed for non-APC mode only
            if nfail < len(config["pyfr"]["dt"]) and mesh_adapter.step > 1:
                # the step is decremented twice
                # to avoid restarting from an incorrect solution
                mesh_adapter.decrement()

                # tstart is lowered
                mesh_adapter.t -= mesh_adapter.K * mesh_adapter.dt_sol
                logger.debug(f"t lowered to: {mesh_adapter.t}")

                # reset gmsh and sol files
                mesh_adapter.gmsh_file = mesh_adapter.dgmsh[mesh_adapter.step]
                mesh_adapter.mesh_file = mesh_adapter.dmsh[mesh_adapter.step]
                mesh_adapter.sol_file = mesh_adapter.dsol[mesh_adapter.step]

                # update dt and pseudo_dt
                mesh_adapter.dt = config["pyfr"]["dt"][nfail]
                logger.info(f"dt: {config['pyfr']['dt'][nfail]}")
                mesh_adapter.pseudo_dt = config['pyfr']['pseudo_dt'][nfail]
                logger.info(f"pseudo_dt: {config['pyfr']['pseudo_dt'][nfail]}")
            else:
                print(f"too many failures ({nfail}), process aborted")
                abort = True
                break
        logger.info(f"adaptation step finished after {time.time() - t0:.2f} sec.")

    adap_duration = time.time() - t0
    logger.info(f"\nadaptation finished in {adap_duration:.2f} sec.")

    # extra pyfr execution with the adapted mesh
    if not abort:
        logger.info("\nextra-computation..")
        # in APC mode the current time must be re-incremented
        if mesh_adapter.mode == "APC":
            mesh_adapter.t += mesh_adapter.K * mesh_adapter.dt_sol
            logger.info(f"increment time counter t to: {mesh_adapter.t}")
        # pyfr final execution
        mesh_adapter.extra_step(mesh_adapter.extra_dt,
                                mesh_adapter.extra_dt_sol)
    logger.info(f"\nextra-computation finished in {time.time() - t0:.2f} sec.")

    # create .series file for paraview
    vtk_list = sum(mesh_adapter.dvtk.values(), [])
    series_file = write_series(mesh_adapter.sol_pattern, vtk_list, os.getcwd())
    logger.info(f"\nproduced series file: {series_file}")

    # save error
    error_df = pd.DataFrame(
        mesh_adapter.derror.items(), columns=["step", "error"]
    )
    error_df.to_csv("error.csv", index=False)
    logger.info("produced error file: error.csv")

    # save QoI history
    if mesh_adapter.stat_crit:
        key_list = [f"{f}_{k}"
                    for f, k_l in zip(mesh_adapter.stat_crit["files"],
                                      mesh_adapter.stat_crit["keys"])
                    for k in k_l]
        logger.debug(f"key_list: {key_list}")
        logger.debug(f"QoI_list: {mesh_adapter.QoI_list}")
        qoi_df = pd.DataFrame(mesh_adapter.QoI_list, columns=key_list)
        qoi_df.to_csv("qoi.csv", index=False)
        logger.info("produced QoI file: qoi.csv")

    # save adaptation steps
    adap_df = pd.DataFrame({"step": mesh_adapter.adaptation_list})
    adap_df.to_csv("adaptation.csv", index=False)
    logger.info("produced adaptation file: adaptation.csv")

    # write execution time statistics
    total_duration = time.time() - t0
    logger.info("\n\n **EXECUTION TIME STATISTICS**")
    logger.info(f"adaptation time: {total_duration:.2f} sec. (100 %)")
    logger.info(f"pyfr execution time: {mesh_adapter.t_run:.2f} sec. "
                f"({mesh_adapter.t_run / total_duration * 100:.2f} %)")
    logger.info(f"pyfr import time: {mesh_adapter.t_import:.2f} sec. "
                f"({mesh_adapter.t_import / total_duration * 100:.2f} %)")
    logger.info(f"pyfr export time: {mesh_adapter.t_export:.2f} sec. "
                f"({mesh_adapter.t_export / total_duration * 100:.2f} %)")
    logger.info(f"pyfr interpolation time: {mesh_adapter.t_interpolate:.2f} sec. "
                f"({mesh_adapter.t_interpolate / total_duration * 100:.2f} %)")
    logger.info(f"vtk extension time: {mesh_adapter.t_extvtk:.2f} sec. "
                f"({mesh_adapter.t_extvtk / total_duration * 100:.2f} %)")
    logger.info(f"madlib time: {mesh_adapter.t_madlib:.2f} sec. "
                f"({mesh_adapter.t_madlib / total_duration * 100:.2f} %)")


if __name__ == '__main__':
    main()
