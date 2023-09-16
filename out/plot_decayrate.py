import matplotlib.pyplot as plt
import pandas as pd

# 读取文件内容
with open("./particles.csv", "r") as f:
    lines = f.readlines()

# 使用 "Time step" 来分隔事件
events = []
event = []
for line in lines:
    if "Time step" in line and event:
        events.append(event)
        event = []
    event.append(line.strip())
if event:  # 添加最后一个事件
    events.append(event)

# 初始化一个DataFrame来存储所有事件的数据
all_data_list = []

for event in events:
    event_data = event[1:]  # 去掉第一行（列名）
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

# 合并所有事件的数据
all_data_df = pd.concat(all_data_list)

# 对相同的时间步进行平均
avg_data = all_data_df.groupby("Time step").mean().reset_index()
final_rate = avg_data["DecayRate"].iloc[-1]
print("rate:=", final_rate)
# 绘制结果
plt.figure(figsize=(10, 6))
plt.plot(
    0.2 * avg_data["Time step"],
    avg_data["DecayRate"],
    label=r"$dN_{\delta}/{N_\Delta}$",
    color="green",
)
plt.plot(
    0.2 * avg_data["Time step"],
    avg_data["AverageRatio"],
    label=r"$\Gamma / \gamma$",
    color="red",
)
# plt.plot(0.02 * avg_data["Time step"], avg_data["PROTON"], label="PROTON", color="blue")
# plt.plot(0.02 * avg_data["Time step"], avg_data["PION"], label="PION", color="green")
# plt.plot(0.02 * avg_data["Time step"], avg_data["DELTA"], label="DELTA", color="red")
plt.xlabel("Time step")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.05, 20)
plt.ylim(0.01, 100)
plt.ylabel("Values")
plt.title("Average DecayRate and Average Ratio vs Time Step")
plt.legend()
plt.grid(True)
plt.savefig("DecayRate.png")
plt.show()
