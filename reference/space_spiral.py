import numpy as np
import matplotlib.pyplot as plt


def space_spiral(a, b, c, theta_max, num_points, kernel_radius, spiral_angle):
    """生成三维空间螺旋线的点。

    参数：
    a, b : 螺旋线在xy平面的公式参数。
    c : 沿z轴的增长速率。
    theta_max : 最大角度（弧度）。
    num_points : 生成的点数。
    kernel_radius : 内核半径。
    spiral_angle : 螺旋角（弧度）。

    返回：
    x, y, z坐标的数组。
    """
    theta = np.linspace(0, theta_max, num_points)
    r = a + b * theta
    x = (r + kernel_radius) * np.cos(theta)
    y = (r + kernel_radius) * np.sin(theta)
    z = c * theta
    return x, y, z


def divide_spiral(x, y, z, num_segments):
    """将螺旋线等分成指定数量的段，并输出各段的节点位移和段号。

    参数：
    x, y, z : 螺旋线的坐标数组。
    num_segments : 要分割成的段数。

    返回：
    每个段的节点位移和段号。
    """
    num_points = len(x)
    segment_length = num_points // num_segments
    segment_displacements = []
    segment_endpoints = []
    for i in range(num_segments):
        start_index = i * segment_length
        end_index = min((i + 1) * segment_length, num_points)
        displacement = np.sqrt(
            (x[end_index - 1] - x[start_index]) ** 2
            + (y[end_index - 1] - y[start_index]) ** 2
            + (z[end_index - 1] - z[start_index]) ** 2
        )
        segment_displacements.append((displacement, i))
        if i < num_segments - 1:
            segment_endpoints.append(
                (x[end_index - 1], y[end_index - 1], z[end_index - 1])
            )
    return segment_displacements, segment_endpoints


# 螺旋线参数
a = 0  # 初始半径
b = 0.5  # 半径增加速率
c = 0.1  # z轴增长速率
theta_max = 4 * np.pi  # 最大角度
num_points = 1000  # 点数
num_segments = 10  # 分成的段数
kernel_radius = 1.0  # 内核半径
spiral_angle = np.pi / 4  # 螺旋角

# 生成螺旋线
x, y, z = space_spiral(a, b, c, theta_max, num_points, kernel_radius, spiral_angle)

# 分割螺旋线
segment_displacements, segment_endpoints = divide_spiral(x, y, z, num_segments)

# 输出各段的节点位移和段号
for displacement, segment_number in segment_displacements:
    print(f"Segment {segment_number}: Displacement = {displacement}")

# 绘图显示螺旋线
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
ax.plot(x, y, z, label="Space Spiral")
ax.scatter(x, y, z, color="red", s=10)  # 显示所有节点
for point in segment_endpoints:
    ax.scatter(*point, color="blue", s=50)  # 显示段的终点
ax.set_title("3D Space Spiral with Points")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.legend()
plt.show()
