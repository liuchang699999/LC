import numpy as np
import math


def read_file():
    # 定义数组和矩阵
    x = np.zeros(50)
    y = np.zeros(50)
    z = np.zeros(50)
    pj = np.zeros(20)
    pf = np.zeros((2, 20))
    je = np.zeros((3, 50), dtype=int)
    jn = np.zeros((6, 50), dtype=int)
    jpj = np.zeros(20, dtype=int)
    jpf = np.zeros((3, 20), dtype=int)
    jeai = np.zeros(20, dtype=int)

    # 读取标题行
    title = input().split()  # 假设标题行输入格式为 20 个连续的四字符字符串

    # 打印标题行
    print("       " + "".join(title))

    # 读取其他输入数据
    nj, njj, n, ne, npj, npf = map(int, input().split())

    # 打印读入的参数
    print(
        f"    NJ= {nj}    NJJ= {njj}    N= {n}    NE= {ne}    NPJ= {npj}    NPF= {npf}"
    )

    print("NO(1) (2) (3) (4) (5) (6)       X         Y         Z")

    # 读取节点信息
    for i in range(nj):
        data = input().split()
        jn[:, i] = list(map(int, data[:6]))
        x[i], y[i], z[i] = map(float, data[6:9])

    # 打印节点信息
    for i in range(nj):
        print(
            f"({i+1:2d}),{jn[0, i]:6d},{jn[1, i]:6d},{jn[2, i]:6d},{jn[3, i]:6d},{jn[4, i]:6d},{jn[5, i]:6d}    {x[i]:10.3f}    {y[i]:10.3f}    {z[i]:10.3f}"
        )

    # 读取节点坐标
    for i in range(nj, nj + njj):
        x[i], y[i], z[i] = map(float, input().split())

    # 打印节点坐标
    for i in range(nj, nj + njj):
        print(f"({i+1:2d}),    {x[i]:10.3f}    {y[i]:10.3f}    {z[i]:10.3f}")

    print("ELEMENT NO.NODE-1 NODE-2 NODE-VMATERIALS")

    # 读取单元信息
    for i in range(ne):
        je[:, i], jeai[i] = map(int, input().split())

    # 打印单元信息
    for i in range(ne):
        print(
            f"{i+1:4d},    {je[0, i]:3d}    {je[1, i]:3d}    {je[2, i]:3d}    {jeai[i]:3d}"
        )

    # 读取节点负载
    if npj != 0:
        print("        NODEL LOADS")
        print(" NO.DISPVALUE")
        for i in range(npj):
            jpj[i], pj[i] = map(float, input().split())
            print(f"        {jpj[i]:7d} {pj[i]:16.3f}")

    # 读取非节点负载
    if npf != 0:
        print("        NON-NODEL LOADS")
        print("NO.ENO.LOAD.MODEL SURFACE")
        for i in range(npf):
            jpf[:, i], pf[:, i] = map(float, input().split())
            print(
                f"{jpf[0, i]:3d}    {jpf[1, i]:3d}    {jpf[2, i]:3d}    {pf[0, i]:10.3f}    {pf[1, i]:10.3f}"
            )


def mke(ke, ie, je, jeai, eai, x, y, z, al, d_d):
    pi = 3.14159
    ii = je[0, ie - 1]
    jj = je[1, ie - 1]
    mt = jeai[ie - 1]

    barL = np.sqrt(
        (x[jj - 1] - x[ii - 1]) ** 2
        + (y[jj - 1] - y[ii - 1]) ** 2
        + (z[jj - 1] - z[ii - 1]) ** 2
    )
    al[ie - 1] = barL

    eai[0, mt - 1] = 2.0e11
    eai[1, mt - 1] = 7.672311e10
    eai[2, mt - 1] = 9.0 * pi * d_d[ie - 1] ** 2 / 25.0
    eai[3, mt - 1] = 369.0 * pi * d_d[ie - 1] ** 4 / 20000.0
    eai[4, mt - 1] = 369.0 * pi * d_d[ie - 1] ** 4 / 40000.0
    eai[5, mt - 1] = 369.0 * pi * d_d[ie - 1] ** 4 / 40000.0

    a1 = eai[0, mt - 1] * eai[2, mt - 1] / barL
    a2 = eai[0, mt - 1] * eai[5, mt - 1] / barL**3
    a3 = eai[0, mt - 1] * eai[5, mt - 1] / barL**2
    a4 = eai[0, mt - 1] * eai[4, mt - 1] / barL**3
    a5 = eai[0, mt - 1] * eai[4, mt - 1] / barL**2
    a6 = eai[1, mt - 1] * eai[3, mt - 1] / barL
    a7 = eai[0, mt - 1] * eai[4, mt - 1] / barL
    a8 = eai[0, mt - 1] * eai[5, mt - 1] / barL

    ke[0, 0] = a1
    ke[0, 6] = -a1
    ke[1, 1] = 12 * a2
    ke[1, 5] = 6 * a3
    ke[1, 7] = -12 * a2
    ke[1, 11] = 6 * a3
    ke[2, 2] = 12 * a4
    ke[2, 4] = -6 * a5
    ke[2, 8] = -12 * a4
    ke[2, 10] = -6 * a5
    ke[3, 3] = a6
    ke[3, 9] = -a6
    ke[4, 4] = 4 * a7
    ke[4, 8] = 6 * a5
    ke[4, 10] = 2 * a7
    ke[5, 5] = 4 * a8
    ke[5, 7] = -6 * a3
    ke[5, 11] = 2 * a8
    ke[6, 6] = a1
    ke[7, 7] = 12 * a2
    ke[7, 11] = -6 * a3
    ke[8, 8] = 12 * a4
    ke[8, 10] = 6 * a5
    ke[9, 9] = a6
    ke[10, 10] = 4 * a7
    ke[11, 11] = 4 * a8

    for i in range(12):
        for k in range(i, 12):
            ke[k, i] = ke[i, k]

    return ke, al


def mr(r, ie, je, x, y, z):
    i = je[0, ie - 1] - 1  # Adjust index for Python's 0-based indexing
    j = je[1, ie - 1] - 1
    k = je[2, ie - 1] - 1

    # Calculate barL and directional cosines
    barL = np.sqrt((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2)
    lxx = (x[j] - x[i]) / barL
    lxy = (y[j] - y[i]) / barL
    lxz = (z[j] - z[i]) / barL

    # Calculate other elements of the directional cosine matrix
    yz = (y[k] - y[i]) * (z[k] - z[j]) - (z[k] - z[i]) * (y[k] - y[j])
    zx = -((x[k] - x[i]) * (z[k] - z[j]) - (z[k] - z[i]) * (x[k] - x[j]))
    xy = (x[k] - x[i]) * (y[k] - y[j]) - (y[k] - y[i]) * (x[k] - x[j])

    l2 = np.sqrt(yz**2 + zx**2 + xy**2)
    lzx = yz / l2
    lzy = zx / l2
    lzz = xy / l2

    # Compute lyx, lyy, lyz
    s1 = (
        (1 - lxx**2) * (x[k] - x[i])
        - lxx * lxy * (y[k] - y[i])
        - lxx * lxz * (z[k] - z[i])
    )
    s2 = (
        (1 - lxy**2) * (y[k] - y[i])
        - lxy * lxx * (x[k] - x[i])
        - lxy * lxz * (z[k] - z[i])
    )
    s3 = (
        (1 - lxz**2) * (z[k] - z[i])
        - lxz * lxx * (x[k] - x[i])
        - lxz * lxy * (y[k] - y[i])
    )

    l3 = np.sqrt(s1**2 + s2**2 + s3**2)
    lyx = s1 / l3
    lyy = s2 / l3
    lyz = s3 / l3

    # Initialize the matrix r with zeros
    r.fill(0)

    # Set values in the transformation matrix
    for ii in range(0, 10, 3):
        r[ii, ii] = lxx
        r[ii, ii + 1] = lxy
        r[ii, ii + 2] = lxz
        r[ii + 1, ii] = lyx
        r[ii + 1, ii + 1] = lyy
        r[ii + 1, ii + 2] = lyz
        r[ii + 2, ii] = lzx
        r[ii + 2, ii + 1] = lzy
        r[ii + 2, ii + 2] = lzz

    return r


########################################
def make(ke, r):
    rt = tran(r)
    tmp = mulv(rt, ke)
    ake = mulv(tmp, r)
    return ake


def tran(r):
    return r.T


def mulv(a, b):
    return np.dot(a, b)


def calm(m, jn, je, ie):
    for i in range(6):
        m[i] = jn[i, je[0, ie] - 1]
        m[i + 6] = jn[i, je[1, ie] - 1]
    return m


def mk(k, ake, m):
    for i in range(12):
        for j in range(12):
            if m[i] != 0 and m[j] != 0:
                k[m[i] - 1, m[j] - 1] += ake[i, j]
    return k


def pe(fe, ip, jpf, pf, al):
    a = pf[0, ip - 1]
    c = pf[1, ip - 1]
    barL = al[jpf[0, ip - 1] - 1]
    ind = jpf[1, ip - 1]
    fe = np.zeros(12)
    if ind == 1:
        fe[1] = (7 * a / 20 + 3 * c / 20) * barL
        fe[5] = (a / 20 + c / 30) * barL**2
        fe[7] = (3 * a / 20 + 7 * c / 20) * barL
        fe[11] = -(a / 30 + c / 20) * barL**2
    # 更多条件的处理需要根据实际情况翻译

    if jpf[2, ip - 1] == 2:
        p = fe[4]
        fe[4] = -fe[5]
        fe[5] = -p
        p = fe[1]
        fe[1] = fe[2]
        fe[2] = p
        p = fe[10]
        fe[10] = -fe[11]
        fe[11] = -p
        p = fe[7]
        fe[7] = fe[8]
        fe[8] = p
    return fe


def mulv12(a, b):
    c = np.zeros(12)
    for i in range(12):
        for j in range(12):
            c[i] += a[i, j] * b[j]
    return c


def mf(p, afe, m):
    for i in range(12):
        if m[i] != 0:
            p[m[i] - 1] += afe[i]
    return p


def slov(ak, p, n):
    d = p.copy()
    for k in range(n - 1):
        for i in range(k + 1, n):
            c = -ak[k, i] / ak[k, k]
            for j in range(i, n):
                ak[i, j] += c * ak[k, j]
            d[i] += c * d[k]
    d[-1] /= ak[-1, -1]
    for i in range(n - 2, -1, -1):
        for j in range(i + 1, n):
            d[i] -= ak[i, j] * d[j]
        d[i] /= ak[i, i]
    return d


def made(ie, jn, je, d):
    ade = np.zeros(12)
    for i in range(6):
        if jn[i, je[0, ie] - 1] != 0:
            ade[i] = d[jn[i, je[0, ie] - 1] - 1]
        if jn[i, je[1, ie] - 1] != 0:
            ade[i + 6] = d[jn[i, je[1, ie] - 1] - 1]
    return ade


# 读取文件
nj, njj, n, ne, nm, npj, npf, jn, x, y, z, je, jeai, jpj, pj, jpf, pf = read_file()

# 初始化
k = np.zeros((100, 100))
ke = np.zeros((12, 12))
ake = np.zeros((12, 12))
x = np.zeros(50)
y = np.zeros(50)
z = np.zeros(50)
al = np.zeros(50)
eai = np.zeros((6, 20))
pj = np.zeros(20)
pf = np.zeros((2, 20))
r = np.zeros((12, 12))
p = np.zeros(50)
d_d = np.zeros(20)
fe = np.zeros(12)
d = np.zeros(50)
ade = np.zeros(12)
de = np.zeros(12)
rt = np.zeros((12, 12))
f = np.zeros(6)
afe = np.zeros(12)
af = np.zeros(12)

je = np.zeros((3, 50), dtype=int)
jn = np.zeros((6, 50), dtype=int)
jpj = np.zeros(20, dtype=int)
jpf = np.zeros((3, 20), dtype=int)
m = np.zeros(12, dtype=int)
jeai = np.zeros(40, dtype=int)


# 给定参数
e = 2.0e11
pi = math.pi
g = 7.672311e10


# 构建总刚度矩阵
for ie in range(1, ne + 1):
    mke(ke, ie, je, jeai, eai, x, y, z, al, d_d)
    mr(r, ie, je, x, y, z)
    make(ke, r, ake)
    calm(m, ie, jn, je)
    mk(k, ake, m)

# 计算节点力
for ip in range(1, npf + 1):
    mr(r, jpf[0][ip], je, x, y, z)
    rt = np.transpose(r)
    pe(fe, ip, jpf, pf, al)
    mulv12(rt, fe, afe)
    calm(m, jpf[0][ip], jn, je)
    mf(p, afe, m)

# 形成节点载荷
for i in range(1, npj + 1):
    p[jpj[i - 1]] += pj[i - 1]

# 求解位移
d = np.linalg.solve(k, p)

# 输出结果
print("RESULTS OF CALCULATION")
print("DISPLACEMENT")
print("NO.E        DX          DY          DZ          RX          RY          RZ")
for kk in range(1, nj + 1):
    f = np.zeros(6)
    for ii in range(1, 7):
        i1 = jn[ii - 1][kk - 1]
        if i1 > 0:
            f[ii - 1] = d[i1 - 1]
    print(
        f"{kk:2d} {f[0]:12.4f} {f[1]:12.4f} {f[2]:12.4f} {f[3]:12.4f} {f[4]:12.4f} {f[5]:12.4f}"
    )

print("单元号       内力DMAX       临界应力D_E       直径D_D")
for ie in range(1, ne + 1):
    made(ie, jn, je, d, ade)
    mulv12(r, ade, de)
    ii = je[0][ie - 1]
    jj = je[1][ie - 1]
    mt = jeai[ie - 1]
    barl = np.sqrt(
        (x[jj - 1] - x[ii - 1]) ** 2
        + (y[jj - 1] - y[ii - 1]) ** 2
        + (z[jj - 1] - z[ii - 1]) ** 2
    )
    de7 = de[6]
    de1 = de[0]
    de2 = de[1]
    de3 = de[2]
    de5 = de[4]
    de6 = de[5]
    de4 = de[3]
    de8 = de[7]
    de12 = de[11]
    de10 = de[9]
    de9 = de[8]
    de11 = de[10]
    d1max = (de7 - de1) * e / barl
    d2max = 3.0 * d_d[ie - 1] * e * (de8 - de2) / barl**2
    tmd11 = d_d[ie - 1] * e * (de12 + 2.0 * de6) / barl
    d2max += tmd11
    d3max = 3.0 * d_d[ie - 1] * e * (de9 - de3) / barl**2
    tmd11 = d_d[ie - 1] * e * (de11 + 2.0 * de5) / barl
    d3max += tmd11
    dwmax = np.sqrt(d2max**2 + d3max**2)
    dtmax = 0.5 * d_d[ie - 1] * g * (de10 - de4) / barl
    dmax = np.sqrt((d1max + dwmax) ** 2 + 3 * dtmax**2)
    d_e = 41 * pi**2 * e * d_d[ie - 1] ** 2 / (1600 * barl**2)
    print(f"{ie:2d} {dmax:12.4f} {d_e:12.4f} {d_d[ie-1]:12.4f}")
