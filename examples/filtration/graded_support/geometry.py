"""Connectivity/coordinate audit of material labels exported by actual OpenLB."""
from collections import deque

import numpy as np

def audit(path,raw,protocol=None):
    rows=np.atleast_1d(np.genfromtxt(path,delimiter=",",names=True,dtype=np.int64))
    points={(int(r["i"]),int(r["j"]),int(r["k"])) for r in rows}
    if len(points)!=len(rows):
        raise RuntimeError("Repeated geometry lattice IDs.")
    dx=raw["dx_m"]; nx=raw["axial_intervals"]
    ny=round(raw["mesh_span_m"]/dx)
    if not points or any(not (0<=i<=nx and 0<=j<ny and 0<=k<ny) for i,j,k in points):
        raise RuntimeError("Invalid exported geometry coordinates.")
    counts=np.bincount(rows["i"],minlength=nx+1)
    if np.any(counts==0) or counts[0]!=raw["inlet_nodes"] or counts[-1]!=raw["outlet_nodes"]:
        raise RuntimeError("Empty pore section or incorrect boundary geometry.")
    for material,i in ((3,0),(4,nx)):
        actual=rows[rows["material"]==material]
        if len(actual)!=raw["inlet_nodes"] or not np.all(actual["i"]==i):
            raise RuntimeError("Pressure boundary material is misplaced.")
    profile_check={}
    if protocol is not None:
        from examples.filtration.graded_support.model import profile
        x,radii=profile(protocol,raw["bulk_radius_m"])
        radii_by_plane=np.interp(np.arange(nx+1)*dx-raw["lead_length_each_m"],x,radii)/dx
        radial=(np.arange(ny)+.5-ny/2)**2
        expected=(radial[None,:,None]+radial[None,None,:]<radii_by_plane[:,None,None]**2).sum(axis=(1,2))
        inside=(rows["j"]+.5-ny/2)**2+(rows["k"]+.5-ny/2)**2<radii_by_plane[rows["i"]]**2
        if not np.array_equal(expected,counts) or not np.all(inside):
            raise RuntimeError("Actual material voxels differ from the declared radius profile.")
        profile_check={"declared_profile_matches_material_labels":True,"profile_planes_checked":nx+1}
    unseen=points.copy(); first=next(iter(unseen)); unseen.remove(first)
    queue=deque([first]); visited=0
    while queue:
        i,j,k=queue.popleft(); visited+=1
        for p in ((i-1,j,k),(i+1,j,k),(i,j-1,k),(i,j+1,k),(i,j,k-1),(i,j,k+1)):
            if p in unseen:
                unseen.remove(p); queue.append(p)
    if unseen:
        raise RuntimeError("Exported fluid geometry contains disconnected components.")
    return {**profile_check,"six_neighbor_connected_components":1,"connected_fluid_nodes":visited,
            "minimum_cross_section_nodes":int(counts.min()),"inlet_nodes":int(counts[0]),
            "outlet_nodes":int(counts[-1]),"coordinate_contract":"zero-based i,j,k; x=i*dx, y=(j+.5)*dx, z=(k+.5)*dx; metres"}
