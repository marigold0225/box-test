import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.special import kn

# Constants
T = 0.155
M_proton = 0.938
M_pion = 0.135
M_delta = 1.232


# Functions to compute rho values
def compute_rho_proton(mu_p):
    y2_proton = kn(2, M_proton / T)
    rho_p = (
        2
        / (2 * np.pi) ** 3
        * 4
        * np.pi
        * T
        * M_proton**2
        * y2_proton
        * np.exp(mu_p / T)
    )
    return rho_p / 0.197**3


def compute_rho_pion(mu_pi):
    y2_pion = kn(2, M_pion / T)
    rho_pi = (
        3 / (2 * np.pi) ** 3 * 4 * np.pi * T * M_pion**2 * y2_pion * np.exp(mu_pi / T)
    )
    return rho_pi / 0.197**3


def compute_rho_delta(mu_delta):
    y2_delta = kn(2, (M_proton + M_pion) / T)
    rho_delta = (
        4.69
        / (2 * np.pi) ** 3
        * 4
        * np.pi
        * T
        * (M_proton + M_pion) ** 2
        * y2_delta
        * np.exp(mu_delta / T)
    )
    return rho_delta / 0.197**3


# Read the data from particles.csv
with open("./particles.csv", "r") as f:
    lines = f.readlines()

events = []
event = []
for line in lines:
    if "Time step" in line and event:
        events.append(event)
        event = []
    event.append(line.strip())
if event:
    events.append(event)

all_data_list = []
for event in events:
    event_data = event[1:]
    df = pd.DataFrame(
        [x.split(",") for x in event_data],
        columns=["Time step", "PROTON", "PION", "DELTA", "DecayRate", "AverageRatio"],
    )
    df = df.astype(
        {
            "Time step": "int32",
            "PROTON": "int32",
            "PION": "int32",
            "DELTA": "int32",
            "DecayRate": "float64",
            "AverageRatio": "float64",
        }
    )
    all_data_list.append(df)

all_data_df = pd.concat(all_data_list)
avg_data = all_data_df.groupby("Time step").mean().reset_index()

# Extract the final values for each particle type
final_proton = avg_data["PROTON"].iloc[-1]
final_pion = avg_data["PION"].iloc[-1]

# Calculate mu values
mu_p_solution = fsolve(lambda mu: compute_rho_proton(mu) * 8000 - final_proton, 0.938)
mu_pi_solution = fsolve(lambda mu: compute_rho_pion(mu) * 8000 - final_pion, 0.135)

# Calculate rho_delta
mu_delta = mu_p_solution + mu_pi_solution
rho_delta_value = compute_rho_delta(mu_delta) * 8000
rho_p = compute_rho_proton(mu_p_solution) * 8000

print("rho_p:=", final_proton, "rho_delta:=", rho_delta_value)
print("mu_p:=", mu_p_solution, "mu_pion:=", mu_pi_solution)
# Plotting
plt.plot(0.2 * avg_data["Time step"], avg_data["PROTON"], label="PROTON", color="blue")
plt.plot(0.2 * avg_data["Time step"], avg_data["PION"], label="PION", color="green")
plt.plot(0.2 * avg_data["Time step"], avg_data["DELTA"], label="DELTA", color="red")
plt.xlabel("Time step")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.05, 20)
plt.ylabel("Average Particle number")
plt.title("Average Particle number vs Time Step")
plt.legend()
plt.grid(True)

# Display the computed rho_delta on the plot
# plt.text(0.05, max(avg_data["DELTA"])/2, f"rho_delta: {rho_delta_value[0]:.2f}", fontsize=10, color="red")
ax = plt.gca()
x_max = ax.get_xlim()[1]
original_ylim = ax.get_ylim()  # 获取原始的纵坐标范围

# 在图的右侧纵坐标对应的位置使用两个加粗的圆点标记 rho_delta_value 和 rho_proton
plt.scatter(
    [x_max] * 2, [rho_delta_value[0], rho_p[0]], color=["red", "blue"], s=100, zorder=5
)

# 在纵坐标刻度线的右侧显示 rho_delta_value 和 rho_proton 的值
ax.text(
    x_max + 0.1,
    rho_delta_value[0],
    f"{rho_delta_value[0]:.2f}",
    verticalalignment="center",
    color="red",
)
ax.text(
    x_max + 0.1, rho_p[0], f"{rho_p[0]:.2f}", verticalalignment="center", color="blue"
)
ratio = rho_delta_value[0] / rho_p[0]

# 在图中的适当位置添加这个比值
plt.text(
    0.1,
    original_ylim[1] - 0.7 * (original_ylim[1] - original_ylim[0]),
    f"N_delta/N_proton = {ratio:.2f}",
    fontsize=10,
)

# 添加新的刻度标签
current_ticks = list(ax.get_yticks())
new_ticks = current_ticks + [rho_delta_value[0], rho_p[0]]
ax.set_yticks(new_ticks)

ax.set_ylim(original_ylim)

plt.savefig("number.png")
plt.show()
