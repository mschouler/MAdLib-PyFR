import gmsh
import json
import numpy as np

from mesh_utils import build_bl, get_mesh_parser, write_mesh


def from_dat(file: str, header_len: int = 2, scale: float = 1) -> list[list[float]]:
    """
    Returns the cleaned list of data points.
    >> dat_fil: path to input_geometry.dat.
    >> pts:  the geometry coordinates in the original referential.
       pts = [[x0, y0, z0], [x1, y1, z1], ..., [xN, yN, zN]]
       where N is the number of points describing the geometry
       and (z0, ..., zN) are null or identical.
    Note:
        When a geometry is closed, the .dat file may contain a redundancy i.e. the last
        point is also the first one in the list. This can become problematic for the
        compressor blade case but has no effect in other cases. Such duplicates are hence removed.
    """
    dat_file = [line.strip() for line in open(file, "r").read().splitlines()]
    pts = [list(map(float, line.split(" "))) for line in dat_file[header_len:]]
    pts = pts[:-1] if pts[0] == pts[-1] else pts
    pts = [[p[0], p[1], 0.] for p in pts] if len(pts[0]) == 2 else pts
    return pts if scale == 1 else [[coord * scale for coord in p] for p in pts]


def reorder_blade(pts: list[list[float]]) -> list[list[float]]:
    """
    **Returns** the blade profile after reordering.
    """
    d = np.sqrt([x**2 + y**2 for x, y, _ in pts])
    start = np.argmin(d)
    if pts[start + 1][1] > pts[start][1]:
        return [[p[0], p[1], p[2]] for p in pts[start:] + pts[:start]]
    else:
        return [[p[0], p[1], p[2]] for p in pts[:start] + pts[start:]]


def build_cylinder_field(
        radius: float, VIn: float, VOut: float,
        XAxis: float, XCenter: float,
        YAxis: float, YCenter: float,
        ZAxis: float = 0.
) -> int:
    """
    **Builds** a cylinder field in the computational domain.
    """
    f_cyl = gmsh.model.mesh.field.add('Cylinder')
    gmsh.model.mesh.field.setNumber(f_cyl, 'Radius', radius)
    gmsh.model.mesh.field.setNumber(f_cyl, 'VIn', VIn)
    gmsh.model.mesh.field.setNumber(f_cyl, 'VOut', VOut)
    gmsh.model.mesh.field.setNumber(f_cyl, 'XAxis', XAxis)
    gmsh.model.mesh.field.setNumber(f_cyl, 'XCenter', XCenter)
    gmsh.model.mesh.field.setNumber(f_cyl, 'YAxis', YAxis)
    gmsh.model.mesh.field.setNumber(f_cyl, 'YCenter', YCenter)
    gmsh.model.mesh.field.setNumber(f_cyl, 'ZAxis', ZAxis)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_cyl)
    return f_cyl


def build_minaniso_field(tag: list[int]):
    """
    **Builds** a MinAniso field in the computational domain.
    """
    f_minaniso = gmsh.model.mesh.field.add('MinAniso')
    gmsh.model.mesh.field.setNumbers(f_minaniso, 'FieldsList', tag)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_minaniso)


def build_2dmesh(pts: list[list[float]], mesh_config: dict, structured: bool) -> list[int]:
    """
    **Builds** the surface mesh of the computational domain
    and returns the surface tag list.
    """
    # domain
    extrusion_layers = mesh_config["domain"].get("extrusion_layers", 0)
    doutlet = mesh_config["domain"].get("outlet", 6.3e-2)
    nodes_inlet = mesh_config["mesh"].get("nodes_inlet", 25)
    nodes_outlet = mesh_config["mesh"].get("nodes_outlet", 17)
    snodes_inlet = mesh_config["mesh"].get("side_nodes_inlet", 31)
    snodes_outlet = mesh_config["mesh"].get("side_nodes_outlet", 31)
    c_snodes = mesh_config["mesh"].get("curved_side_nodes", 7)
    # wall
    le = mesh_config["mesh"].get("le", 16)
    te = mesh_config["mesh"].get("te", 16)
    nodes_sp2 = mesh_config["mesh"].get("nodes_sp2", 42)
    nodes_sp3 = mesh_config["mesh"].get("nodes_sp3", 42)
    nodes_sp4 = mesh_config["mesh"].get("nodes_sp4", 14)
    nodes_sp7 = mesh_config["mesh"].get("nodes_sp7", 57)
    nodes_sp8 = mesh_config["mesh"].get("nodes_sp8", 32)
    # boundary layer
    bl = mesh_config["boundary_layer"].get("bl", False)
    size = mesh_config["boundary_layer"].get("size", 5e-5)
    ratio = mesh_config["boundary_layer"].get("ratio", 1.15)
    thickness = mesh_config["boundary_layer"].get("thickness", 4e-3)
    sizefar = mesh_config["boundary_layer"].get("sizefar", 5e-4)
    # wake refinement
    wake = mesh_config["wake"].get("wake", False)
    cyl_vin = mesh_config["wake"].get("cyl_vin", 8e-4)
    cyl_vout = mesh_config["wake"].get("cyl_vout", 5e-3)
    cyl_xaxis = mesh_config["wake"].get("cyl_xaxis", 1.675e-2)
    cyl_xcenter = mesh_config["wake"].get("cyl_xcenter", 8.364e-2)

    wall = reorder_blade(pts)
    pt_wall = [gmsh.model.geo.addPoint(p[0], p[1], p[2]) for p in wall]

    # blade splines and transfinite curves
    spl_1 = gmsh.model.geo.addSpline(pt_wall[:35])
    spl_2 = gmsh.model.geo.addSpline(pt_wall[35 - 1:88])
    spl_3 = gmsh.model.geo.addSpline(pt_wall[88 - 1:129])
    spl_4 = gmsh.model.geo.addSpline(pt_wall[129 - 1:157])
    spl_5 = gmsh.model.geo.addSpline(pt_wall[157 - 1:168])
    spl_6 = gmsh.model.geo.addSpline(pt_wall[168 - 1:179])
    spl_7 = gmsh.model.geo.addSpline(pt_wall[179 - 1:245])
    spl_8 = gmsh.model.geo.addSpline(pt_wall[245 - 1:287])
    spl_9 = gmsh.model.geo.addSpline(pt_wall[287 - 1:322] + [pt_wall[0]])
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_1, le // 2, "Progression", 1.02)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_2, nodes_sp2, "Progression", 1.03)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_3, nodes_sp3, "Progression", 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_4, nodes_sp4, "Progression", 0.94)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_5, te // 2, "Progression", 0.97)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_6, te // 2, "Progression", 1.025)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_7, nodes_sp7, "Progression", 1.015)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_8, nodes_sp8, "Progression", 0.955)
    gmsh.model.geo.mesh.setTransfiniteCurve(spl_9, le // 2, "Progression", 0.9)
    spl_list = [spl_1, spl_2, spl_3, spl_4, spl_5, spl_6, spl_7, spl_8, spl_9]
    blade_loop = gmsh.model.geo.addCurveLoop(spl_list)

    # domain construction points
    pt_323 = gmsh.model.geo.addPoint(-6e-2, -5e-2, 0.)
    pt_324 = gmsh.model.geo.addPoint(-1.5e-2, -2.5e-2, 0.)
    pt_325 = gmsh.model.geo.addPoint(0., -1.7e-2, 0.)
    pt_326 = gmsh.model.geo.addPoint(1.270973e-02, -1.164466e-02, 0.)
    pt_327 = gmsh.model.geo.addPoint(2.585445e-02, -7.360298e-03, 0.)
    pt_328 = gmsh.model.geo.addPoint(3.934429e-02, -4.053609e-03, 0.)
    pt_329 = gmsh.model.geo.addPoint(5.308943e-02, -1.631280e-03, 0.)
    pt_330 = gmsh.model.geo.addPoint(6.7e-2, 0., 0.)
    pt_331 = gmsh.model.geo.addPoint(6.7e-2 + doutlet, 0., 0.)
    pt_332 = gmsh.model.geo.addPoint(6.7e-2 + doutlet, 4.039e-2, 0.)
    pt_333 = gmsh.model.geo.addPoint(6.7e-2, 4.039e-2, 0.)
    pt_334 = gmsh.model.geo.addPoint(5.308943e-02, 3.875872e-02, 0.)
    pt_335 = gmsh.model.geo.addPoint(3.934429e-02, 3.633639e-02, 0.)
    pt_336 = gmsh.model.geo.addPoint(2.585445e-02, 3.302970e-02, 0.)
    pt_337 = gmsh.model.geo.addPoint(1.270973e-02, 2.874534e-02, 0.)
    pt_338 = gmsh.model.geo.addPoint(0., 2.339e-2, 0.)
    pt_339 = gmsh.model.geo.addPoint(-1.5e-2, 1.539e-2, 0.)
    pt_340 = gmsh.model.geo.addPoint(-6e-2, -9.61e-3, 0.)

    # domain construction lines
    l_10 = gmsh.model.geo.addLine(pt_340, pt_323, tag=10)
    l_11 = gmsh.model.geo.addLine(pt_323, pt_324, tag=11)
    l_12 = gmsh.model.geo.addLine(pt_339, pt_340, tag=12)
    l_13 = gmsh.model.geo.addLine(pt_324, pt_339, tag=13)
    l_14 = gmsh.model.geo.addLine(pt_324, pt_325, tag=14)
    l_15 = gmsh.model.geo.addLine(pt_325, pt_326, tag=15)
    l_16 = gmsh.model.geo.addLine(pt_326, pt_327, tag=16)
    l_17 = gmsh.model.geo.addLine(pt_327, pt_328, tag=17)
    l_18 = gmsh.model.geo.addLine(pt_328, pt_329, tag=18)
    l_19 = gmsh.model.geo.addLine(pt_329, pt_330, tag=19)
    l_20 = gmsh.model.geo.addLine(pt_330, pt_331, tag=20)
    l_21 = gmsh.model.geo.addLine(pt_331, pt_332, tag=21)
    l_22 = gmsh.model.geo.addLine(pt_332, pt_333, tag=22)
    l_23 = gmsh.model.geo.addLine(pt_333, pt_334, tag=23)
    l_24 = gmsh.model.geo.addLine(pt_334, pt_335, tag=24)
    l_25 = gmsh.model.geo.addLine(pt_335, pt_336, tag=25)
    l_26 = gmsh.model.geo.addLine(pt_336, pt_337, tag=26)
    l_27 = gmsh.model.geo.addLine(pt_337, pt_338, tag=27)
    l_28 = gmsh.model.geo.addLine(pt_338, pt_339, tag=28)

    # transfinite curves on non-blade boundaries
    gmsh.model.geo.mesh.setTransfiniteCurve(l_10, nodes_inlet, "Progression", 1.)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_13, nodes_inlet, "Progression", 1.)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_21, nodes_outlet, "Progression", 1.)
    # bottom / top periodicity
    bottom_tags = [l_11, l_14, l_15, l_16, l_17, l_18, l_19, l_20]
    top_tags = [l_12, l_28, l_27, l_26, l_25, l_24, l_23, l_22]
    # bottom non-curved side nodes
    gmsh.model.geo.mesh.setTransfiniteCurve(l_11, snodes_inlet, "Progression", 1.)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_20, snodes_outlet, "Progression", 1.)
    # bottom curved side nodes
    _ = [gmsh.model.geo.mesh.setTransfiniteCurve(l_i, c_snodes, "Progression", 1.)
         for l_i in bottom_tags[1:-1]]
    # periodic boundaries /y direction
    gmsh.model.geo.synchronize()
    translation = [1, 0, 0, 0, 0, 1, 0, 0.04039, 0, 0, 1, 0, 0, 0, 0, 1]
    for tid, bid in zip(top_tags, bottom_tags):
        gmsh.model.mesh.setPeriodic(1, [tid], [bid], translation)

    # closed curve loop and computational domain surface definition
    cloop_2 = gmsh.model.geo.addCurveLoop([l_10, l_11, l_12, l_13])
    cloop_3 = gmsh.model.geo.addCurveLoop(
        [-l_13, l_14, l_15, l_16, l_17, l_18, l_19, l_20,
            l_21, l_22, l_23, l_24, l_25, l_26, l_27, l_28,]
    )
    surf_1 = gmsh.model.geo.addPlaneSurface([cloop_2], tag=1002)
    surf_2 = gmsh.model.geo.addPlaneSurface([cloop_3, blade_loop], tag=1003)
    gmsh.model.geo.mesh.setTransfiniteSurface(surf_1)

    # Fields definition
    # Boundary Layer
    if bl:
        build_bl(spl_list, size, ratio, thickness, sizefar, structured)
    # Cylinder #1
    if wake:
        f_cyl1 = build_cylinder_field(
            9e-3,
            cyl_vin,
            cyl_vout,
            cyl_xaxis,
            cyl_xcenter,
            -1.171e-3,
            1.9754e-2
        )
        # Cylinder #2
        f_cyl2 = build_cylinder_field(
            1.62e-2, 1.6e-3, 5e-3, 2.01e-2, 8.699e-2, -1.406e-3, 1.9519e-2
        )
        # MinAniso
        build_minaniso_field([f_cyl1, f_cyl2])

    # define physical groups for boundary conditions
    surf_tag = [surf_1, surf_2]
    if extrusion_layers == 0:
        gmsh.model.geo.addPhysicalGroup(2, surf_tag, tag=100, name="fluid")
        gmsh.model.geo.addPhysicalGroup(1, [l_10], tag=10, name="inlet")
        print(f"2D BC: Inlet tags are {[l_10]}")
        gmsh.model.geo.addPhysicalGroup(1, [l_21], tag=20, name="outlet")
        print(f"2D BC: Outlet tags are {[l_21]}")
        gmsh.model.geo.addPhysicalGroup(1, spl_list, tag=30, name="wall")
        print(f"2D BC: Wall tags are {spl_list}")
        gmsh.model.geo.addPhysicalGroup(1, top_tags, tag=40, name="periodic_0_l")
        print(f"2D BC: Top tags are {top_tags}")
        gmsh.model.geo.addPhysicalGroup(1, bottom_tags, tag=50, name="periodic_0_r")
        print(f"2D BC: Bottom tags are {bottom_tags}")
    else:
        gmsh.model.geo.addPhysicalGroup(2, surf_tag, tag=100, name="periodic_1_l")
    return surf_tag


def build_3dmesh(mesh_config: dict, surf_tag: list[int]):
    """
    **Performs** an extrusion along the z axis.
    - h_size (float): the total extruded depth.
    """
    h_size = mesh_config["domain"]["extrusion_size"]
    extrusion_layers = mesh_config["domain"]["extrusion_layers"]
    print(f"DEBUG -- surf_tag: {surf_tag}")
    ext_tag = [gmsh.model.geo.extrude(
        [(2, s)], 0, 0, h_size, [extrusion_layers], [1], True) for s in surf_tag]
    # retrieve extruded surfaces and volumes
    vol = [tu[-1] for tu in [ext_tag[0][1], ext_tag[1][1]]]
    top = [tu[-1] for tu in [ext_tag[0][0], ext_tag[1][0]]]
    # 1st block
    inlet = [ext_tag[0][2][-1]]
    perlo = [ext_tag[0][3][-1]]
    perup = [ext_tag[0][5][-1]]
    # 2nd block
    perlo += [tu[-1] for tu in ext_tag[1][3:10]]
    outlet = [ext_tag[1][10][-1]]
    perup += [tu[-1] for tu in ext_tag[1][11:18]]
    wall = [tu[-1] for tu in ext_tag[1][18:]]
    # create physical groups
    gmsh.model.geo.addPhysicalGroup(3, vol, tag=2000, name="fluid")
    print("3D BC: vol tag is 2000")
    gmsh.model.geo.addPhysicalGroup(2, top, tag=2001, name="periodic_1_r")
    print("3D BC: periodic_span_l tag is 100")
    print("3D BC: periodic_span_r tag is 2001")
    gmsh.model.geo.addPhysicalGroup(2, perlo, tag=2003, name="periodic_0_l")
    print("3D BC: periodic_vert_l tag is 2003")
    gmsh.model.geo.addPhysicalGroup(2, perup, tag=2004, name="periodic_1_r")
    print("3D BC: periodic_vert_r tag is 2004")
    gmsh.model.geo.addPhysicalGroup(2, inlet, tag=2005, name="inlet")
    print("3D BC: inlet tag is 2005")
    gmsh.model.geo.addPhysicalGroup(2, outlet, tag=2006, name="outlet")
    print("3D BC: outlet tag is 2006")
    gmsh.model.geo.addPhysicalGroup(2, wall, tag=2007, name="wall")
    print("3D BC: wall tag is 2007")


def main():
    """
    This program automatically generates
    a simple cylinder in a rectangle mesh with gmsh.

    Note: the domain dimensions are hardcoded.
    """
    parser = get_mesh_parser()
    parser.add_argument(
        "-in", "--input", type=str,
        required=True, help="path to the geometry"
    )
    args = parser.parse_args()

    # all mesh parameters
    if args.config:
        with open(args.config) as jfile:
            mesh_config = json.load(jfile)["mesh_config"]
    else:
        mesh_config = {
            "domain": {
                "extrusion_layers": 0,
                "extrusion_size": 1.4e-2,
                "outlet": 6.3e-2
            },
            "mesh": {},
            "boundary_layer": {},
            "wake": {}
        }

    # initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.logger.start()
    gmsh.model.add("model")

    # load surface
    pts = from_dat(args.input)
    surf_tag = build_2dmesh(pts, mesh_config, args.structured)
    if args.structured:
        gmsh.model.geo.mesh.setRecombine(2, abs(surf_tag))
    gmsh.model.geo.synchronize()
    if mesh_config["domain"]["extrusion_layers"] > 0:
        _3d = True
        build_3dmesh(mesh_config, surf_tag)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
    else:
        _3d = False
        gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(args.order)

    try:
        gmsh.fltk.run()
    except Exception:
        print("WARNING -- graphical rendering failed")

    write_mesh(args.name, args.format, args.outdir, args.log, args.geo, _3d)


if __name__ == "__main__":
    main()
