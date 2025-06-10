import argparse
import gmsh
import os


def build_bl(
        blade_tag: list[int],
        size: float,
        ratio: float,
        thickness: float,
        sizefar: float = -1,
        structured: bool = False
):
    """
    Builds the boundary layer:

    - wall_tag: the tag of the boundary layer surface
    - size: the size of the first layer
    - ratio: the growth ratio from one layer to the next one
    - thickness: the overall thickness of the boundary layer mesh
    - structured: whether it should be made of quads or tris
    """
    f_bl = gmsh.model.mesh.field.add('BoundaryLayer')
    gmsh.model.mesh.field.setNumbers(f_bl, 'CurvesList', blade_tag)
    gmsh.model.mesh.field.setNumber(f_bl, 'Size', size)
    gmsh.model.mesh.field.setNumber(f_bl, 'Ratio', ratio)
    gmsh.model.mesh.field.setNumber(f_bl, 'Thickness', thickness)
    if sizefar > 0:
        gmsh.model.mesh.field.setNumber(f_bl, 'SizeFar', sizefar)
    gmsh.model.mesh.field.setNumber(f_bl, 'Quads', int(structured))
    gmsh.model.mesh.field.setAsBoundaryLayer(f_bl)


def get_mesh_kwd(mesh: list[str], kwd: str) -> int:
    try:
        idx = next(idx for idx, el in enumerate(mesh) if kwd in el)
    except StopIteration:
        raise Exception(f"ERROR -- no '{kwd}' entry in mesh file")
    return idx


def get_mesh_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--config", type=str, default="", help="mesh config file")
    parser.add_argument("-n", "--name", type=str, default="madlib_pyfr_mesh", help="mesh name")
    parser.add_argument(
        "-out", "--outdir", type=str, default=os.getcwd(), help="mesh output directory")
    parser.add_argument("-f", "--format", type=str, default="msh", help="mesh format")
    parser.add_argument("-k", "--order", type=int, default=1, help="mesh order")
    parser.add_argument("-l", "--log", action="store_true", help="generate gmsh .log file")
    parser.add_argument("-g", "--geo", action="store_true", help="generate gmsh .geo file")
    parser.add_argument(
        "-s", "--structured", action="store_true", help="generates a structured mesh")
    return parser


def reformat_2d(mesh: list[str]) -> list[str]:
    """
    Fixes gmsh default .mesh format in 2D.
    """
    idx = get_mesh_kwd(mesh, "Dimension")
    mesh[idx] = " Dimension 2"
    del mesh[idx + 1]

    vert_idx = get_mesh_kwd(mesh, "Vertices")
    n_vert = int(mesh[vert_idx + 1])
    for id in range(vert_idx + 2, vert_idx + 2 + n_vert):
        line_data = list(map(float, mesh[id].split()))
        mesh[id] = " " * 4 + f"{line_data[0]:>20}" + \
                   " " * 4 + f"{line_data[1]:>20}" + \
                   " " * 4 + f"{int(line_data[-1]):>20}"
    return mesh


def write_mesh(name: str, format: str, outdir: str, log: bool, geo: bool, _3d: bool = False):
    """
    Writes the mesh name.format to outdir.
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        print(f"created {outdir} repository")

    if geo:
        print(f"writing {name}.geo_unrolled to {outdir}")
        gmsh.write(os.path.join(outdir, name + ".geo_unrolled"))

    print(f"writing {name}.{format} to {outdir}")
    gmsh.write(os.path.join(outdir, name + f".{format}"))

    # medit formatting
    if format == "mesh":
        print(f"medit formatting of {name}.mesh")
        mesh_file = os.path.join(outdir, name + ".mesh")
        mesh = open(mesh_file, "r").read().splitlines()
        mesh = reformat_2d(mesh) if not _3d else mesh
        with open(mesh_file, 'w') as ftw:
            ftw.write("\n".join(mesh))

    if log:
        log_content = gmsh.logger.get()
        log_file = open(os.path.join(outdir, name + ".log"), "w")
        print(f"writing {name}.log to {outdir}")
        log_file.write("\n".join(log_content))
    gmsh.logger.stop()
    gmsh.finalize()
