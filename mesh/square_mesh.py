import gmsh
import json

from mesh_utils import get_mesh_parser, write_mesh


def build_2dmesh(mesh_config: dict, structured: bool) -> int:
    """
    Builds the 2d mesh and returns the tag of the computational domain.

    - mesh_config: a dictionary containing domain, mesh and boundary layer parameters
    - structured: whether to build a structured (True) or unstructured (False) mesh
    """
    # extract domain parameters
    center = mesh_config["domain"]["center"]
    width = mesh_config["domain"]["width"]
    height = mesh_config["domain"]["height"]

    # domain corners
    bot_left = gmsh.model.geo.addPoint(center[0] - width, center[1] - height, 0.)
    bot_right = gmsh.model.geo.addPoint(center[0] + width, center[1] - height, 0.)
    top_left = gmsh.model.geo.addPoint(center[0] - width, center[1] + height, 0.)
    top_right = gmsh.model.geo.addPoint(center[0] + width, center[1] + height, 0.)

    # domain boundaries
    left = gmsh.model.geo.addLine(top_left, bot_left)
    bottom = gmsh.model.geo.addLine(bot_left, bot_right)
    right = gmsh.model.geo.addLine(bot_right, top_right)
    top = gmsh.model.geo.addLine(top_right, top_left)

    # extract mesh parameters
    inodes = mesh_config["mesh"]["inlet_nodes"]
    snodes = mesh_config["mesh"]["side_nodes"]

    # parameterized boundaries
    gmsh.model.geo.mesh.setTransfiniteCurve(left, inodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(right, inodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, snodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, snodes, "Progression", 1)

    # periodic boundaries /y direction
    gmsh.model.geo.synchronize()
    translation = [1, 0, 0, 0, 0, 1, 0, 20, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [top], [bottom], translation)

    # periodic boundaries /x direction
    gmsh.model.geo.synchronize()
    translation = [1, 0, 0, 20, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [right], [left], translation)

    # closed curve loop and computational domain surface definition
    boundary_loop = gmsh.model.geo.addCurveLoop([left, bottom, right, top])
    surf_tag = gmsh.model.geo.addPlaneSurface([boundary_loop], tag=1000)

    gmsh.model.geo.mesh.setTransfiniteSurface(surf_tag)

    # define physical groups for boundary conditions
    gmsh.model.geo.addPhysicalGroup(1, [right], tag=2, name="periodic_0_r")
    # gmsh.model.geo.addPhysicalGroup(1, [right], tag=2, name="rightside")
    print(f"BC: Right tags are {[right]}")
    gmsh.model.geo.addPhysicalGroup(1, [left], tag=3, name="periodic_0_l")
    # gmsh.model.geo.addPhysicalGroup(1, [left], tag=3, name="leftside")
    print(f"BC: Left tags are {[left]}")
    gmsh.model.geo.addPhysicalGroup(1, [top], tag=4, name="periodic_1_r")
    print(f"BC: Top tags are {[top]}")
    gmsh.model.geo.addPhysicalGroup(1, [bottom], tag=5, name="periodic_1_l")
    print(f"BC: Bottom tags are {[bottom]}")
    gmsh.model.geo.addPhysicalGroup(2, [surf_tag], tag=1, name="Fluid")
    return surf_tag


def main():
    """
    This program automatically generates a simple cylinder in a rectangle mesh with gmsh.

    Note: the domain dimensions are hardcoded.
    """
    parser = get_mesh_parser()
    args = parser.parse_args()

    # all mesh parameters
    if args.config:
        with open(args.config) as jfile:
            mesh_config = json.load(jfile)
    else:
        mesh_config = {
            "domain": {
                "center": (0, 0),
                "width": 10,
                "height": 10
            },
            "mesh": {
                "inlet_nodes": 21,
                "side_nodes": 21
            },
        }

    # initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.logger.start()
    gmsh.model.add("model")

    surf_tag = build_2dmesh(mesh_config, args.structured)
    if args.structured:
        gmsh.model.geo.mesh.setRecombine(2, abs(surf_tag))
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.model.mesh.setOrder(args.order)

    try:
        gmsh.fltk.run()
    except Exception:
        print("WARNING -- graphical rendering failed")

    write_mesh(args.name, args.format, args.outdir, args.log, args.geo)


if __name__ == "__main__":
    main()
