import json
import logging
import os
import time

from subprocess import run
from typing import Callable, Any, TypeVar

T = TypeVar('T')


def set_logger(log_name: str, log_level: int):
    """
    Sets the logger.
    """
    logger = logging.getLogger()
    logger.setLevel(log_level)
    ft = logging.Formatter(
        '%(asctime)s : %(name)s : %(levelname)s : %(message)s'
    )
    # file handler
    file_handler = logging.FileHandler(log_name, mode="w")
    file_handler.setFormatter(ft)
    logger.addHandler(file_handler)
    # console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(ft)
    logger.addHandler(console_handler)


def run_cmd(
        exec_cmd: list[str], output_file: str, tout: float | None = None
):
    """
    Wrapper around the subprocess.run function.
    """
    with open(os.path.join(output_file), "wb") as out:
        run(exec_cmd, stdout=out, stderr=out, check=True, timeout=tout)


def get_files(outdir: str, ext: str, next: str = "") -> list[str]:
    """
    Retuns the timestamp ordered list of files in outdir
    that contain (do not contain) the ext (next) pattern in their name.
    """
    list_dir = [os.path.join(outdir, f) for f in os.listdir(outdir)
                if os.path.isfile(os.path.join(outdir, f)) and ext in f]
    if next:
        list_dir = [f for f in list_dir if next not in f]
    file_format = list_dir[-1].split(".")[-1]
    sorted_list = sorted(
        list_dir,
        key=lambda x: float(x.split('_')[-1][:-(len(file_format) + 1)])
    )
    return sorted_list


def write_series(sol_pattern: str, sol_list: list[str], outdir: str) -> str:
    """
    Write a .series file and returns its path.
    """
    content: dict = {"file-series-version": "1.0."}
    content["files"] = [
        {"name": sname, "time": sid} for sid, sname in enumerate(sol_list)
    ]
    json_content = json.dumps(content, indent=4)
    outfile = os.path.join(outdir, f"{sol_pattern}.vtu.series")
    with open(outfile, 'w') as ftw:
        ftw.write(json_content)
    return os.path.join(outdir, f"{sol_pattern}.vtu.series")


def time_function(f: Callable[..., T]) -> Callable[..., tuple[T, float]]:  # type: ignore  # noqa
    """
    Timing function wrapper.
    """
    def wrapper(self, *args: Any, **kwargs: Any) -> tuple[T, float]:
        start_time = time.time()
        result = f(self, *args, **kwargs)
        return result, time.time() - start_time
    return wrapper
