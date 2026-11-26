"""Standalone scientific figure for the coordinator's temperature-support test."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent / "cryo-window-extension"
df = pd.read_csv(root / "windows.csv")
selected = df.query("scenario=='high' and sar=='clamped'")
supports = [(-150, 0), (-123, 0), (-123, -35), (-100, -40)]
labels = ["Entire trajectory", "−123 to 0 °C", "−123 to −35 °C", "−100 to −40 °C"]
fig, ax = plt.subplots(figsize=(9, 4.3), layout="constrained")
y = np.arange(4)
for design, offset, color, label in [("margin", -.18, "#21618C", "Published ramp"),
                                     ("equal", .18, "#BA4A00", "Constant, equal integrated duty")]:
    values = [selected.query("design==@design and lower==@lo and upper==@hi").rate.iloc[0]
              for lo, hi in supports]
    ax.barh(y+offset, values, height=.32, color=color, label=label)
    for pos, value in zip(y+offset, values):
        ax.text(value+.5, pos, f"{value:.2f}", va="center", fontsize=9)
ax.axvline(50, color="#444444", ls="--", lw=1)
ax.set(yticks=y, yticklabels=labels, xlim=(0, 84), xlabel="Minimum modeled warming rate (°C/min)")
ax.invert_yaxis()
ax.set_title("Changing the rate-screen interval reverses the schedule ranking", loc="left", fontsize=13)
ax.legend(loc="upper left", bbox_to_anchor=(0, -.18), ncol=2, frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.text(.20, -.13, "Vitrify phantom · 0.50 mm · 1.83 mg Fe/mL muscle · clamped SAR · adiabatic\nThermal diagnostics only; ice, mechanics and survival are not modeled.", fontsize=9)
for suffix in ["png", "svg"]:
    fig.savefig(root / f"criterion_comparison.{suffix}", dpi=180, bbox_inches="tight")
print(root / "criterion_comparison.png")
