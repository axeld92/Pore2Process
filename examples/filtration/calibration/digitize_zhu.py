"""Reproduce Figure 6 M1 readings from the JPEG or PDF image, without fitting.

Four touching diamond pairs are split using visually recorded centres. Other
centres are bounding-box midpoints of the dark connected components. Pixel
uncertainty is retained separately from unreported experimental uncertainty.
"""
from pore2process_demo.paths import project_path
import argparse
from pathlib import Path

import numpy as np
from PIL import Image

from examples.filtration.calibration.run import ROOT, read_json, sha256, write_json


def digitize(config):
    source=project_path(ROOT/config["image"])
    if sha256(source)!=config["image_sha256"]:
        raise ValueError("Figure 6 image changed; repeat the visual digitization audit.")
    a=np.asarray(Image.open(source).convert("RGB"))
    mask=np.max(a,axis=2)<config["dark_threshold"]
    x0,y0,x1,y1=config["pixel_region_xyxy"]
    region=np.zeros_like(mask)
    region[y0:y1,x0:x1]=mask[y0:y1,x0:x1]
    centres=[]
    components=[]
    for y,x in zip(*np.where(region)):
        if not region[y,x]:
            continue
        stack=[(int(y),int(x))]
        region[y,x]=False
        xs,ys=[],[]
        while stack:
            yy,xx=stack.pop()
            xs.append(xx); ys.append(yy)
            for dy,dx in ((-1,0),(1,0),(0,-1),(0,1)):
                yn,xn=yy+dy,xx+dx
                if 0<=yn<region.shape[0] and 0<=xn<region.shape[1] and region[yn,xn]:
                    region[yn,xn]=False; stack.append((yn,xn))
        if len(xs)<=config["minimum_component_pixels"]:
            continue
        box=[min(xs),min(ys),max(xs),max(ys)]
        components.append({"bbox_xyxy":box,"dark_pixels":len(xs)})
        split=next((s for s in config["touching_pairs"] if s["bbox_xyxy"]==box),None)
        if split:
            centres.extend(split["visually_read_centres_xy"])
        else:
            if len(xs)>config.get("maximum_single_component_pixels",90):
                raise ValueError(f"Unreviewed merged component: {box}")
            centres.append([(box[0]+box[2])/2,(box[1]+box[3])/2])
    centres.sort()
    if len(centres)!=config["expected_point_count"]:
        raise ValueError("Marker count differs from the visually audited series.")
    points=[]
    for i,(x,y) in enumerate(centres,1):
        Js,Jm=pixels_to_flux(config,x,y)
        points.append({"point_id":i,"pixel_x":x,"pixel_y":y,
                       "support_flux_m_h":Js,"total_flux_LMH":Jm,
                       "experimental_uncertainty":None})
    return {"config":config,"components":sorted(components,key=lambda c:c["bbox_xyxy"]),"points":points}


def pixels_to_flux(config,x,y,anchors=None):
    a=anchors or config["axis_anchors"]
    Js=(x-a["x0_px"])*100/(a["x100_px"]-a["x0_px"])
    Jm=600+(a["y600_px"]-y)*3000/(a["y600_px"]-a["y3600_px"])
    return Js,Jm


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--source",choices=("jpeg","pdf"),default="jpeg")
    args=parser.parse_args()
    config=read_json(ROOT/'literature/sources/zhu2021/extraction.json')["digitization_pdf" if args.source=="pdf" else "digitization"]
    write_json(args.output,digitize(config))
    print(f"Saved {config['expected_point_count']} M1 readings to {args.output}")


if __name__=="__main__":
    main()
