# -*- coding: utf-8 -*-
"""
论文2 图表生成脚本 fig_paper2.py
Fig 2  分解瀑布（1990-2023，3因子+4因子）
Fig 3  覆盖情景柱状图（年可规避 DALY）
Fig 4  连续剂量反应曲线 + 政策达标区（数据：gbd_step3c_doseresp_floored.csv）
Fig 5  2024-2050 投影轨迹（数据：gbd_step4_bapc_traj.csv）

路径按 WSL 实际位置修改 DATA_DIR 即可。全部英文标注（WSL matplotlib 无中文字体问题）。
运行：python3 fig_paper2.py   → 输出 fig2.pdf fig3.pdf fig4.pdf fig5.pdf 到当前目录
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.bbox": "tight",
})

DATA_DIR = "/mnt/f/HRS/GBD/work"   # <-- 按需修改
OUT = "."

# ============================================================
# Fig 2. Decomposition waterfall（锁版数字，来自 gbd_step2b_decomp3）
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4))

def waterfall(ax, start, steps, labels, unit="million DALYs", title="", total_lab="",
              hi_index=None, hi_color="#e08214"):
    """steps: list of (label, value). start: baseline value.
    hi_index: 0-based index of the step bar to highlight (e.g., the PAF-change bar)."""
    cum = start
    xs, bottoms, heights, colors = [], [], [], []
    xs.append(0); bottoms.append(0); heights.append(start); colors.append("#bdbdbd")
    for i, (lab, v) in enumerate(steps):
        xs.append(i + 1)
        if v >= 0:
            bottoms.append(cum); heights.append(v)
            colors.append(hi_color if i == hi_index else "#2166ac")
        else:
            bottoms.append(cum + v); heights.append(-v); colors.append("#b2182b")
        cum += v
    xs.append(len(steps) + 1)
    bottoms.append(0); heights.append(cum); colors.append("#636363")
    ax.bar(xs, heights, bottom=bottoms, color=colors, width=0.6, edgecolor="white")
    for x, b, h in zip(xs, bottoms, heights):
        ax.text(x, b + h + max(heights) * 0.02, f"{h:.2f}", ha="center",
                va="bottom", fontsize=7.5)
    # 汇总标注放在终点柱正上方
    ax.text(xs[-1], heights[-1] + max(heights) * 0.09, total_lab,
            ha="center", fontsize=9, fontweight="bold", color="#333333")
    ax.set_xticks(xs)
    ax.set_xticklabels(["1990"] + labels + ["2023"], fontsize=7.5)
    ax.set_ylabel(unit, fontsize=8.5)
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, heights[-1] * 1.18)

# 左：全因痴呆 3 因子（百万）
ax = axes[0]
waterfall(ax, 13.46,
          [("Population", 8.85), ("Ageing", 13.12), ("Rate", 4.28)],
          ["Pop.\ngrowth", "Ageing\n(structure)", "DALY\nrate"],
          title="All-cause dementia DALYs",
          total_lab="+26.25\nmillion")

# 右：抑郁归因 4 因子（百万）；PAF change 柱(第4步, index 3)暖色高亮
ax = axes[1]
waterfall(ax, 0.085,
          [("Population", 0.068), ("Ageing", 0.090), ("Rate", 0.030), ("PAF", 0.060)],
          ["Pop.\ngrowth", "Ageing\n(structure)", "DALY\nrate", "PAF\nchange"],
          title="Depression-attributable dementia DALYs",
          total_lab="+0.248\nmillion",
          hi_index=3)

axes[1].text(0.99, 0.03, "Note different y-axis scales\nbetween panels",
             transform=axes[1].transAxes, ha="right", va="bottom", fontsize=6.5,
             color="#555555")

fig.tight_layout()
fig.savefig(f"{OUT}/fig2_decomp.pdf"); fig.savefig(f"{OUT}/fig2_decomp.png")
plt.close(fig)

# ============================================================
# Fig 4 (合并原 Fig3+Fig4). 剂量反应曲线 + 情景柱状图, 上下双面板, 共用 coverage x 轴
# ============================================================
df = pd.read_csv(f"{DATA_DIR}/gbd_step3c_doseresp_floored.csv")
# v2: 三 η 线全画（数据 csv 已含 0.2/0.3/0.5）
g_all = df[df["panel"] == "Global"]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.6),
                               gridspec_kw={"height_ratios": [3, 2]})

# --- Panel A: 连续剂量反应曲线（三 η 线：0.2 / 0.3 primary / 0.5）---
eta_style = {0.2: ("#999999", 1.3, "--"), 0.3: ("#2166ac", 2.2, "-"), 0.5: ("#b2182b", 1.3, "--")}
for eta in [0.2, 0.3, 0.5]:
    ge = g_all[g_all["eta"] == eta].sort_values("c_target")
    col, lw, ls = eta_style[eta]
    lab = "eta = 0.3 (primary)" if eta == 0.3 else f"eta = {eta}"
    ax1.plot(ge["c_target"] * 100, ge["avoid_daly"] / 1e3,
             color=col, lw=lw, ls=ls, label=lab)
ax1.axvspan(13.65, 50, color="#f4a582", alpha=0.3, label="Policy target zone")
ax1.axvline(13.65, color="#d6604d", ls="--", lw=1)
ax1.text(13.65, ax1.get_ylim()[0], " 1.5× WHO 2030", rotation=90, va="bottom",
         fontsize=7, color="#d6604d")
ax1.axvline(50, color="#d6604d", ls="--", lw=1)
ax1.text(50, ax1.get_ylim()[0], " S1 = 50%", rotation=90, va="bottom",
         fontsize=7, color="#d6604d")
c_now_global = 9.1   # mdd_mat 全球 2021 实际值, 投稿前核对
ax1.axvline(c_now_global, color="#636363", ls=":", lw=1.2)
ax1.text(c_now_global + 1, ax1.get_ylim()[1] * 0.05, "current global\ncoverage (9.1%)",
         fontsize=7, color="#636363")
ax1.axhline(0, color="black", lw=0.5)
ax1.set_ylabel("Thousands of DALYs averted per year", fontsize=8.5)
ax1.set_title("A  Dose–response: coverage vs avoidable dementia burden", fontsize=9)
ax1.legend(frameon=False, fontsize=8, loc="upper left")
ax1.tick_params(labelsize=8)
ax1.set_xlim(0, 100)

# --- Panel B: 三个政策情景柱状图 (与 A 共用 x 轴含义) ---
scen_x = [13.65, 50, 100]
scen_lab = ["1.5× WHO 2030\n(13.65%)", "50% coverage\n(S1)", "Universal\ncoverage"]
# v2: 从 gbd_step3c_doseresp_floored.csv 动态读数（不再硬编码；锁版 S1=35,690）
_targets = [0.1365, 0.50, 1.00]
_g3 = df[(df["eta"] == 0.3) & (df["panel"] == "Global")].sort_values("c_target")
# v2.1: np.interp 沿曲线插值（13.65% 不在网格上，argmin 会抓到相邻网格行）
val = [float(np.interp(t, _g3["c_target"].values, _g3["avoid_daly"].values))
       for t in _targets]
cols = ["#92c5de", "#4393c3", "#2166ac"]
bars = ax2.bar(scen_x, [v / 1e3 for v in val], width=11, color=cols, edgecolor="white")
for b, v in zip(bars, val):
    ax2.text(b.get_x() + b.get_width() / 2, v / 1e3 + 18, f"{v / 1e3:,.0f}",
             ha="center", fontsize=8)
ax2.set_xlim(0, 100)
ax2.set_ylim(0, max([v / 1e3 for v in val]) * 1.15)
ax2.set_xticks(scen_x)
ax2.set_xticklabels(scen_lab, fontsize=8)
ax2.set_ylabel("Thousands per year (η = 0.3)", fontsize=8.5)
ax2.set_xlabel("Treatment coverage (%)", fontsize=8.5)
ax2.set_title("B  Annual avoidable burden at policy scenario coverage levels", fontsize=9)
ax2.tick_params(axis="y", labelsize=8)

fig.tight_layout()
fig.savefig(f"{OUT}/fig34_doseresp_scenarios.pdf")
fig.savefig(f"{OUT}/fig34_doseresp_scenarios.png")
plt.close(fig)

# ============================================================
# Fig 5. 2024-2050 投影轨迹（数据驱动）
# ============================================================
traj = pd.read_csv(f"{DATA_DIR}/gbd_step4_bapc_traj.csv")  # year, attrib, avoidB, avoidC
traj = traj.sort_values("year")

fig, ax = plt.subplots(figsize=(4.8, 3.4))
ax.plot(traj["year"], traj["attrib"] / 1e4, color="#636363", lw=2,
        label="Attributable DALYs (status quo)")
ax.fill_between(traj["year"], traj["avoidB"] / 1e4, color="#4393c3", alpha=0.25,
                label="Avoided: 50% coverage by 2035")
ax.fill_between(traj["year"], traj["avoidC"] / 1e4, color="#2166ac", alpha=0.25,
                label="Avoided: universal coverage by 2035")
ax.plot(traj["year"], (traj["attrib"] - traj["avoidB"]) / 1e4, color="#4393c3", lw=1.5)
ax.plot(traj["year"], (traj["attrib"] - traj["avoidC"]) / 1e4, color="#2166ac", lw=1.5)
ax.set_xlabel("Year")
ax.set_ylabel("Thousands of DALYs per year")
ax.set_title("Projected depression-attributable dementia DALYs, 2024–2050")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/fig5_projection.pdf"); fig.savefig(f"{OUT}/fig5_projection.png")
plt.close(fig)

print("done: fig2-fig5 (.pdf + .png)")
