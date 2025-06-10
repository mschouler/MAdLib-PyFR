import gmsh
import json

from mesh_utils import build_bl, get_mesh_parser, write_mesh


def build_2dmesh(mesh_config: dict, structured: bool) -> int:
    """
    Builds the 2d mesh and returns the tag of the computational domain.

    - mesh_config: a dictionary containing domain,
      mesh and boundary layer parameters
    - structured: whether to build a structured (True)
      or unstructured (False) mesh
    """
    # extract domain parameters
    center = mesh_config["domain"]["center"]
    diam = mesh_config["domain"]["diameter"]
    dinlet = mesh_config["domain"]["dinlet"]
    doutlet = mesh_config["domain"]["doutlet"]
    dbottom = mesh_config["domain"]["dbottom"]
    dtop = mesh_config["domain"]["dtop"]

    # cylinder points
    ccenter = gmsh.model.geo.addPoint(*center, 0.)
    c1 = gmsh.model.geo.addPoint(center[0] - diam / 2, center[1], 0.)
    c2 = gmsh.model.geo.addPoint(center[0], center[1] + diam / 2, 0.)
    c3 = gmsh.model.geo.addPoint(center[0] + diam / 2, center[1], 0.)
    c4 = gmsh.model.geo.addPoint(center[0], center[1] - diam / 2, 0.)

    # domain corners
    bot_left = gmsh.model.geo.addPoint(
        center[0] - dinlet, center[1] - dbottom, 0.
    )
    bot_right = gmsh.model.geo.addPoint(
        center[0] + doutlet, center[1] - dbottom, 0.
    )
    top_left = gmsh.model.geo.addPoint(
        center[0] - dinlet, center[1] + dtop, 0.
    )
    top_right = gmsh.model.geo.addPoint(
        center[0] + doutlet, center[1] + dtop, 0.
    )

    # cylinder boundaries
    ca1 = gmsh.model.geo.addCircleArc(c1, ccenter, c2)
    ca2 = gmsh.model.geo.addCircleArc(c2, ccenter, c3)
    ca3 = gmsh.model.geo.addCircleArc(c3, ccenter, c4)
    ca4 = gmsh.model.geo.addCircleArc(c4, ccenter, c1)

    # domain boundaries
    inlet = gmsh.model.geo.addLine(top_left, bot_left)
    bottom = gmsh.model.geo.addLine(bot_left, bot_right)
    outlet = gmsh.model.geo.addLine(bot_right, top_right)
    top = gmsh.model.geo.addLine(top_right, top_left)

    # extract mesh parameters
    inodes = mesh_config["mesh"]["inlet_nodes"]
    snodes = mesh_config["mesh"]["side_nodes"]
    wnodes = mesh_config["mesh"]["wall_nodes"]

    # parameterized boundaries
    gmsh.model.geo.mesh.setTransfiniteCurve(ca1, wnodes // 4, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(ca2, wnodes // 4, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(ca3, wnodes // 4, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(ca4, wnodes // 4, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(inlet, inodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(outlet, inodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(bottom, snodes, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(top, snodes, "Progression", 1)

    # boundary_layer
    bl = mesh_config["boundary_layer"]["bl"]
    if bl:
        bl_size = mesh_config["boundary_layer"]["bl_size"]
        bl_ratio = mesh_config["boundary_layer"]["bl_ratio"]
        bl_thickness = mesh_config["boundary_layer"]["bl_thickness"]
        build_bl(
            [ca1, ca2, ca3, ca4], bl_size, bl_ratio, bl_thickness, structured=structured
        )

    # closed curve loop and computational domain surface definition
    boundary_loop = gmsh.model.geo.addCurveLoop([inlet, bottom, outlet, top])
    cylinder_loop = gmsh.model.geo.addCurveLoop([ca1, ca2, ca3, ca4])
    surf_tag = gmsh.model.geo.addPlaneSurface(
        [boundary_loop, cylinder_loop], tag=1000
    )

    # define physical groups for boundary conditions
    gmsh.model.geo.addPhysicalGroup(
        1, [ca1, ca2, ca3, ca4], tag=1, name="wall"
    )
    print(f"BC: Wall tags are {[ca1, ca2, ca3, ca4]}")
    gmsh.model.geo.addPhysicalGroup(
        1, [inlet, top, bottom], tag=2, name="inlet"
    )
    print(f"BC: Inlet tags are {[inlet, top, bottom]}")
    gmsh.model.geo.addPhysicalGroup(
        1, [outlet], tag=3, name="outlet"
    )
    print(f"BC: Outlet tags are {[outlet]}")
    gmsh.model.geo.addPhysicalGroup(2, [surf_tag], tag=4, name="Fluid")
    return surf_tag


def main():
    """
    This program automatically generates
    a simple cylinder in a rectangle mesh with gmsh.

    Note: the domain dimensions are hardcoded.
    """
    parser = get_mesh_parser()
    parser.add_argument("-d", "--diameter", type=float, help="cylinder diameter", default=1)
    parser.add_argument("--coarse", action="store_true", help="generate a coarse mesh")
    args = parser.parse_args()

    # all mesh parameters
    if args.config:
        with open(args.config) as jfile:
            mesh_config = json.load(jfile)
    else:
        D = args.diameter
        mesh_config = {
            "domain": {
                "center": (0, 0),
                "diameter": D,
                "dinlet": 8 * D,
                "doutlet": 35 * D,
                "dbottom": 8 * D,
                "dtop": 8 * D,
                "center": (0, 0)
            }
        }
        if args.coarse:
            mesh_config["mesh"] = {
                "inlet_nodes": 15,
                "side_nodes": 20,
                "wall_nodes": 30
            }
            mesh_config["boundary_layer"] = {
                "bl": True,
                "bl_size": 0.2,
                "bl_ratio": 1.1,
                "bl_thickness": 2
            }
        else:
            mesh_config["mesh"] = {
                "inlet_nodes": 30,
                "side_nodes": 80,
                "wall_nodes": 60
            }
            mesh_config["boundary_layer"] = {
                "bl": True,
                "bl_size": 0.1,
                "bl_ratio": 1.1,
                "bl_thickness": 2
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
