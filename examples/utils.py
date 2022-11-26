import numpy as np
import pandas as pd

from spotpy.objectivefunctions import kge, rsquared, bias

cm1 = 1 / 2.54  # centimeters in inches
figw_1c = 8.5 * cm1  # maximum width for 1 column
figw_2c = 17.5 * cm1  # maximum width for 2 columns


def scatter_1(ax, x, y, label="treatment", xlabel="obs", ylabel="sim",
              best_fit=True, veg_ws=None):
    compare = pd.DataFrame({"x": x, "y": y})
    if veg_ws is not None:
        compare[veg_ws == 0] = np.nan
    compare = compare.dropna()
    ax.plot(compare["x"], compare["y"], marker="o",
            linestyle="None", markersize=2, color="k", fillstyle="none")
    ax.plot([-0.1, 10], [-0.1, 10], color="dodgerblue", alpha=0.7,
            linewidth="0.8")
    ax.axes.set_xticks(np.arange(0, 10 + 2, 2))
    ax.axes.set_yticks(np.arange(0, 10 + 2, 2))
    ax.set_xlim(-0.1, 10)
    ax.set_ylim(-0.1, 10)
    if best_fit:
        p = np.polyfit(compare["x"], compare["y"], 1)
        f = np.poly1d(p)

        # Calculating new x's and y's
        x_new = np.linspace(0, 10, y.size)
        y_new = f(x_new)

        # Plotting the best fit line with the equation as a legend in latex
        ax.plot(x_new, y_new, "r--", linewidth="0.8")
    ax.text(0.02, 0.9, f"{label}", color="k", zorder=10,
            transform=ax.transAxes)
    ax.text(0.6, 0.04, "$Bias$ = " + str(
        round(bias(np.asarray(compare["y"]), np.asarray(compare["x"])), 2)) +
            "\n" + "$R^2$ = " + str(
        round(rsquared(np.asarray(compare["y"]), np.asarray(compare["x"])),
              2)) +
            "\n" + "KGE = " + str(
        round(kge(np.asarray(compare["y"]), np.asarray(compare["x"])), 2)),
            color="k", zorder=10, transform=ax.transAxes)
    return ax
