import numpy as np
from OCC.Core.gp import gp_Pnt
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCC.Core.Geom import Geom_BSplineCurve
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Core.BRepTools import BRepTools_WireExplorer
from OCC.Core.TopExp import TopExp_Explorer


class STEP_EXPORTER():
    def __init__(self):
        pass
        """
        self.r = r
        self.z = z
        self.phi = phi
        self.inner_surface_index = inner_surface_index
        self.outer_surface_index = outer_surface_index if outer_surface_index is not None else r.shape[0] - 1

        if self.inner_surface_index is not None:
            if self.inner_surface_index == self.outer_surface_index: # If equal, ignore inner
                print(f"Provided surface indices are equal, proceeding with one outer surface at index {self.outer_surface_index}")
                self.inner_surface_index = None
            if self.inner_surface_index > self.outer_surface_index: # If inner is larger than outer, switch them
                print(f"Provided index for inner surface is larger than outer surface, indices will be swapped")
                self.inner_surface_index, self.outer_surface_index = (self.outer_surface_index, self.inner_surface_index)
"""
    def is_wire_closed(self, wire):
        exp = BRepTools_WireExplorer(wire)
        vertices = []
        while exp.More():
            vertex = exp.CurrentVertex()
            vertices.append(vertex)
            exp.Next()
        return vertices[0].IsSame(vertices[-1])


    def check_edge_face_connectivity(self, solid: TopoDS_Shape):
        """
        Checks the number of faces each edge is connected to.
        Useful for ensuring watertightness.
        """
        edge_face_map = {}

        face_exp = TopExp_Explorer(solid, TopAbs_FACE)
        while face_exp.More():
            face = face_exp.Current()
            edge_exp = TopExp_Explorer(face, TopAbs_EDGE)
            while edge_exp.More():
                edge = edge_exp.Current()

                found = False
                for existing_edge in edge_face_map:
                    if edge.IsSame(existing_edge):
                        edge_face_map[existing_edge] += 1
                        found = True
                        break
                if not found:
                    edge_face_map[edge] = 1

                edge_exp.Next()
            face_exp.Next()

        # Summary
        bad_edges = 0
        for edge, count in edge_face_map.items():
            if count != 2:
                print(f"!! Edge used in {count} face(s)")
                bad_edges += 1

        return bad_edges


    def cylindrical_to_cartesian_scaled(self, r_slice, z_slice, phi_angle, scale=1000):
        x = r_slice * np.cos(phi_angle) * scale
        y = r_slice * np.sin(phi_angle) * scale
        z = z_slice * scale

        return x, y, z


    def create_uniform_bspline_wire(self, x, y, z, degree=3):
        n = len(x)

        if not (np.isclose(x[0], x[-1]) and np.isclose(y[0], y[-1]) and np.isclose(z[0], z[-1])):
            x = np.append(x, x[0])
            y = np.append(y, y[0])
            z = np.append(z, z[0])
            n += 1

        poles = TColgp_Array1OfPnt(1, n)
        for i in range(n):
            poles.SetValue(i + 1, gp_Pnt(x[i], y[i], z[i]))

        num_knots = n - degree + 1
        knots = TColStd_Array1OfReal(1, num_knots)
        mults = TColStd_Array1OfInteger(1, num_knots)
        for i in range(1, num_knots + 1):
            knots.SetValue(i, float(i - 1))
            mults.SetValue(i, 1)
        mults.SetValue(1, degree + 1)
        mults.SetValue(num_knots, degree + 1)

        bspline = Geom_BSplineCurve(poles, knots, mults, degree, False)
        edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
        wire = BRepBuilderAPI_MakeWire(edge).Wire()
        return wire


    def create_wire_list(self, r, z, phi):
        """
        Creates a list of wires at each phi cross-section at a certain surface index.
        """
        wires = []
        num_phi = phi.shape[0]

        for i in range(num_phi):
            try:
                r_slice = r[:, i]
                z_slice = z[:, i]
                phi_angle = float(phi[i])

                xx, yy, zz = self.cylindrical_to_cartesian_scaled(r_slice, z_slice, phi_angle)
                wire = self.create_uniform_bspline_wire(xx, yy, zz)
                closed = self.is_wire_closed(wire)

                if not closed:
                    print(f"Warning: Wire {i} is NOT closed!")

                wires.append(wire)

            except Exception as e:
                print(f"Skipping section {i} due to error: {e}")

        return wires


    def generate_plasma_solid(self, r, z, phi, r_in=None, z_in=None, phi_in=None):
        """
        Builds the solid by lofting through spline wires at each phi cross-section.
        """
        surface_maker = BRepOffsetAPI_ThruSections(True, True, 1e-6)  # solid=True, ruled=True
        num_sec = phi.shape[0]
        added_sections = 0
        print(f"Generating solid...")
        wires = self.create_wire_list(r, z, phi)

        for wire in wires:
            surface_maker.AddWire(wire)
            added_sections += 1

        if all(v is not None for v in (r_in, z_in, phi_in)):
            num_sec *= 2
            wires_in = self.create_wire_list(r_in, z_in, phi_in)
            for wire in wires_in:
                surface_maker.AddWire(wire)
                added_sections += 1

        print(f"Added {added_sections} wire sections out of {num_sec}.")

        surface_maker.Build()
        solid = surface_maker.Shape()

        explorer = TopologyExplorer(solid)
        num_faces = sum(1 for _ in explorer.faces())
        num_edges = sum(1 for _ in explorer.edges())
        num_shells = sum(1 for _ in explorer.shells())
        num_solids = sum(1 for _ in explorer.solids())
        print(f"Solid details: {num_faces} faces, {num_edges} edges, {num_shells} shell(s), {num_solids} solid(s)")
        return solid


    def generate_coil_solid(self, xx, yy, zz):
        scale = 1000  # m to mm
        npts = len(xx[1])
        # Append all section wires
        section_pts_list = []
        for k in range(npts - 1):
            pts = [gp_Pnt(scale * float(xx[p, k]), scale * float(yy[p, k]), scale * float(zz[p, k])) for p in range(4)]
            section_pts_list.append(pts)

        # Append a copy of the first section for closure
        section_pts_list.append(section_pts_list[0])

        # Now build wires
        section_wires = []
        for pts in section_pts_list:
            edges = []
            for m in range(4):
                edge = BRepBuilderAPI_MakeEdge(pts[m], pts[(m + 1) % 4]).Edge()
                edges.append(edge)
                wire_builder = BRepBuilderAPI_MakeWire()
            for edge in edges:
                wire_builder.Add(edge)
            section_wires.append(wire_builder.Wire())

        # Loft solid through cross-section wires
        loft = BRepOffsetAPI_ThruSections(True, True, 1.0e-6)
        for wire in section_wires:
            loft.AddWire(wire)
        loft.Build()
        solid = loft.Shape()
        return solid


    def compound_solid(self, operation, solid=None, compound=None):
        from OCC.Core.TopoDS import TopoDS_Compound
        from OCC.Core.BRep import BRep_Builder

        if operation == "create":
            compound = TopoDS_Compound()
            builder = BRep_Builder()
            builder.MakeCompound(compound)
            return compound
        elif operation == "add":
            builder = BRep_Builder()
            builder.Add(compound, solid)
            return compound
        else:
            raise ValueError("Invalid operation. Use 'create' or 'add'.")



    def check_watertightness(self, solid: TopoDS_Shape):
        analyzer = BRepCheck_Analyzer(solid)
        if not analyzer.IsValid():
            print("!! Solid is not topologically valid.")
            return False

        if solid.ShapeType() != TopAbs_SOLID:
            print("!! Shape is not a TopAbs_SOLID.")
            return False

        explorer = TopologyExplorer(solid)
        shells = list(explorer.shells())

        if len(shells) != 1:
            print("!! Solid contains multiple or no shells.")
            return False

        bad_edges = self.check_edge_face_connectivity(solid)
        if bad_edges != 0:
            print(f"!! Found {bad_edges} edges not shared by two faces.")
            return False

        return True

    def export_step_file(self, solid, filename="output.step"):
      #  if not self.check_watertightness(solid):
       #     print("!! Warning: solid may not be watertight. Proceeding anyway...")

        print(f"Exporting solid to STEP file: {filename}")

        step_writer = STEPControl_Writer()
        step_writer.Transfer(solid, STEPControl_AsIs)
        status = step_writer.Write(filename)

        if status != IFSelect_RetDone:
            raise RuntimeError("STEP export failed.")

        print("STEP export successful.")
