import pandas as pd
import matplotlib.pyplot as plt

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

# 初始化每个时间步的粒子数目为0的DataFrame
all_data = pd.DataFrame(columns=["Time step", "PROTON", "PION", "DELTA"])

# 遍历每个事件并累积粒子数目
for event in events:
    event_data = event[1:]  # 去掉第一行（列名）
    df = pd.DataFrame(
        [x.split(",") for x in event_data],
        columns=["Time step", "PROTON", "PION", "DELTA", "Total Energy"],
    )
    df = df.astype(
        {"Time step": "int32", "PROTON": "int32", "PION": "int32", "DELTA": "int32"}
    )

    # 合并数据，对相同的时间步进行求和
    all_data = (
        pd.concat([all_data, df])
        .groupby("Time step")
        .sum(numeric_only=True)
        .reset_index()
    )

# 计算平均值
avg_data = all_data.copy()
avg_data["PROTON"] /= len(events)
avg_data["PION"] /= len(events)
avg_data["DELTA"] /= len(events)

# 绘制结果
plt.figure(figsize=(10, 6))
plt.plot(0.02 * avg_data["Time step"], avg_data["PROTON"], label="PROTON", color="blue")
plt.plot(0.02 * avg_data["Time step"], avg_data["PION"], label="PION", color="green")
plt.plot(0.02 * avg_data["Time step"], avg_data["DELTA"], label="DELTA", color="red")
plt.xlabel("Time step")
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.05, 20)
plt.ylabel("Average Number of Particles")
plt.title("Average Number of Particles vs Time Step")
plt.legend()
plt.grid(True)
plt.savefig("box_test.png")
plt.show()
