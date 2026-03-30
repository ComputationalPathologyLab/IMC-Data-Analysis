from pathlib import Path
import sys
import tifffile as tiff
import numpy as np

# --------------------------------------------------
# Usage:
#   python stack_one_roi.py ROI001_D13 data/img/ROI001_D13.tiff
#
# If no arguments are provided, it defaults to:
#   ROI001_D13 -> data/img/ROI001_D13.tiff
# --------------------------------------------------

if len(sys.argv) >= 3:
    roi_dir = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
else:
    roi_dir = Path("ROI001_D13")
    out_path = Path("data/img/ROI001_D13.tiff")

out_path.parent.mkdir(parents=True, exist_ok=True)

channel_order = [
    "115In_CD20.ome.tiff",
    "127I_127I.ome.tiff",
    "131Xe_131Xe.ome.tiff",
    "134Xe_134Xe.ome.tiff",
    "138Ba_138Ba.ome.tiff",
    "141Pr_CD68.ome.tiff",
    "142Nd_CD138.ome.tiff",
    "143Nd_CD47.ome.tiff",
    "144Nd_CD31.ome.tiff",
    "145Nd_Tbet.ome.tiff",
    "146Nd_BCL2.ome.tiff",
    "147Sm_CD44.ome.tiff",
    "148Nd_CD163.ome.tiff",
    "149Sm_CD45RO.ome.tiff",
    "150Nd_PDL1.ome.tiff",
    "151Eu_MHC_II.ome.tiff",
    "152Sm_CD66b.ome.tiff",
    "153Eu_LAG3.ome.tiff",
    "154Sm_TIM3.ome.tiff",
    "155Gd_FoxP3.ome.tiff",
    "156Gd_CD4.ome.tiff",
    "158Gd_CTLA4.ome.tiff",
    "159Tb_Brachiury.ome.tiff",
    "160Gd_CD11c.ome.tiff",
    "161Dy_CD86.ome.tiff",
    "162Dy_CD8.ome.tiff",
    "163Dy_ERG.ome.tiff",
    "164Dy_S100.ome.tiff",
    "165Ho_PD1.ome.tiff",
    "167Er_GranzymeB.ome.tiff",
    "169Tm_CD56.ome.tiff",
    "170Er_CD3.ome.tiff",
    "171Yb_CD21.ome.tiff",
    "172Yb_CD10.ome.tiff",
    "173Yb_Ki67.ome.tiff",
    "174Yb_CD72a.ome.tiff",
    "175Lu_GATA3.ome.tiff",
    "176Yb_CD80.ome.tiff",
    "191Ir_DNA1.ome.tiff",
    "193Ir_DNA2.ome.tiff",
    "208Pb_208Pb.ome.tiff",
    "209Bi_aSMA.ome.tiff",
    "80ArAr_80ArAr.ome.tiff",
]

if not roi_dir.exists():
    raise FileNotFoundError(f"ROI directory does not exist: {roi_dir}")

imgs = []
shapes = []

print(f"Reading ROI from: {roi_dir}")
print(f"Writing stacked TIFF to: {out_path}")

for fn in channel_order:
    p = roi_dir / fn
    if not p.exists():
        raise FileNotFoundError(f"Missing channel: {p}")
    img = tiff.imread(p)
    if img.ndim != 2:
        raise ValueError(f"{p} is not 2D. Shape: {img.shape}")
    imgs.append(img)
    shapes.append(img.shape)

if len(set(shapes)) != 1:
    raise ValueError(f"Channel shapes do not match: {set(shapes)}")

stack = np.stack(imgs, axis=0)  # CYX

tiff.imwrite(out_path, stack)

print(f"Wrote {out_path}")
print(f"Shape: {stack.shape}")
print(f"Dtype: {stack.dtype}")
print(f"Channels stacked: {len(channel_order)}")