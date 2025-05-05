import numpy as np

class dot:
    def __init__(self, coords, parco = np.array([0,0,0])):
        x0,x1,x2 = coords[0], coords[1],coords[2]
        self.x = np.array([x0, x1, x2])
        self.coords = coords
        self.my_triengles = []
        self.direction = np.array([0,0,0])
        self.chrl = 0
        self.moovable = True
        self.parentco = parco
    def smt(self):
        if len(self.my_triengles)>0:
            sum_of_s = 0
            for t in self.my_triengles:
                sum_of_s += t.s()
            return sum_of_s
        else:
            return -1

    def eat(self, dots, r, dfr):
        eaten = []
        for d in dots:
            if 0 < np.linalg.norm(d.coords - self.coords) < r and (d not in dfr):
                eaten.append(d)
                for t in d.my_triengles:
                    if t.d1 == d:
                        t.d1 = self
                    if t.d2 == d:
                        t.d2 = self
                    if t.d3 == d:
                        t.d3 = self
                    self.my_triengles.append(t)
                    #print(1)

        coo = np.array([0,0,0])
        pcoo = np.array([0,0,0])
        for d in eaten:
            coo = coo + d.coords / (len(eaten) + 1)
            pcoo = pcoo + d.parentco / (len(eaten) + 1)
        self.coords = coo + self.coords / (len(eaten) + 1)
        self.parentco = pcoo + self.parentco / (len(eaten) + 1)
        return(eaten)


class triengle:
    def __init__(self, d1, d2, d3):
        self.d1 = d1
        self.d2 = d2
        self.d3 = d3
        d1.my_triengles.append(self)
        d2.my_triengles.append(self)
        d3.my_triengles.append(self)

    def s(self):
        l1 = np.dot(self.d1.coords - self.d2.coords, self.d1.coords - self.d2.coords)**0.5
        l2 = np.dot(self.d1.coords - self.d3.coords, self.d1.coords - self.d3.coords) ** 0.5
        l3 = np.dot(self.d2.coords - self.d3.coords, self.d2.coords - self.d3.coords) ** 0.5
        p = (l1 + l2 + l3) / 2
        return (abs(p * (p - l1) * (p - l2) * (p - l3)))**0.5


def proj(a,b):
    return b * np.dot(a,b) / np.linalg.norm(b)**2


def lsht(contur,t):
    return (np.linalg.norm((contur(t) - contur(t + 0.0001)))) / 0.0001

def contur2(t):
    return np.array([1 * np.cos(t), 1 * np.sin(t), 0.5 * np.sin(4 * t)])

def contur(t):
    return np.array([1 * np.sin(t), 1 * np.sin(2 * t), 1 * np.sin(3 * t)])

def contur1(t):
    __ = 4
    k0 = 1/2
    if t <= np.pi / 2:
        return np.array([t, 0, k0 * np.sin(__ * t)])
    if np.pi / 2<=t < np.pi:
        return np.array([np.pi / 2, t - np.pi / 2, k0 * np.sin(__ * t)])
    if np.pi< t <= 3 * np.pi / 2:
        return np.array([np.pi / 2 - (t - np.pi), np.pi / 2, k0 * np.sin(__ * t)])
    if np.pi* 3 / 2 <=t <= np.pi * 2:
        return np.array([0, np.pi / 2 - (t - 3/2 *np.pi), k0 * np.sin(__ * t)])

def razbivka(contur, N):
    ans = []
    tl = 0
    lns = 0
    moment = 0
    dt = 0.001 / N
    for _ in np.linspace(0,2 * np.pi, 10000):
        tl += lsht(contur, _) * 2 * np.pi * 0.0001
    for _ in np.linspace(0,tl,N):
        while lns < _:
            moment += dt
            lns += lsht(contur, moment) * dt
        ans.append(moment)
    return(np.array(ans))


n = 401
T_max = 2 * np.pi - 2 * np.pi / n
pararam = 0.6
L = 30
M = 10
MD = 30
layers = []
triengles = []
layer0 = []
LMAX = 0
for t in np.linspace(0, T_max, n):
    layer0.append(dot(contur(t)))
    if t != 0:
        LMAX += np.linalg.norm(contur(t) - contur(t - T_max / n)) / (n-1)
layers.append(layer0)





def add_layer1(N, prev_layer, lmax, layer_num):
    par = 0.9
    layer = []
    lsr = 0
    nd = None
    nd0 = None

    for i in range(0, N):
        usd = N
        angle = np.pi
        sch = 0

        c0 = prev_layer[i].coords
        while angle > np.pi * par and sch < 10:
            sch += 1
            v1 = prev_layer[(i - sch)%usd].coords - c0
            v2 = prev_layer[(i + sch) % usd].coords - c0
            angle = abs(np.arccos(np.dot(v1,v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))

        vspdot = (prev_layer[(i - sch)%usd].coords + prev_layer[(i + sch) % usd].coords) / 2
        dir = (vspdot - c0) / np.linalg.norm(vspdot - c0)
        newdot = dir * min(lmax, 5 * (np.linalg.norm(prev_layer[(i - 1)%usd].coords - c0) + np.linalg.norm(prev_layer[(i + 1)%usd].coords - c0)) / 2)
        #newdot = dir * lmax
        newdot = c0 + newdot - proj(newdot, prev_layer[(i - sch)%usd].coords - prev_layer[(i + sch) % usd].coords)
        chekvec = c0 - layers[len(layers) - 2][i].coords
        if np.dot(newdot - c0,chekvec) < 0:
            newdot = newdot - 1.5 * (newdot - c0)
        lsr += (np.linalg.norm(prev_layer[(i - 1)%usd].coords - c0) + np.linalg.norm(prev_layer[(i + 1)%usd].coords - c0)) / 2 / N
        pd = nd
        nd = dot(newdot, c0)
        if i == 0:
            nd0 = nd
        layer.append(nd)
        if i != 0:
            triengles.append(triengle(pd, nd, prev_layer[(i - 1)%usd]))
            triengles.append(triengle(prev_layer[i], nd, prev_layer[(i - 1) % N]))




    triengles.append(triengle(nd, nd0, prev_layer[usd-1]))
    triengles.append(triengle(prev_layer[0], nd0, prev_layer[usd - 1]))

    lsr = 0
    dots_for_remove = []
    mid_dots = []
    prev_dot_wasnt_del = True
    pprev_dot_wasnt_del = True
    for d_num in range(len(layer)):
        lsr += np.linalg.norm(layer[d_num].coords - layer[(d_num + 1) % N].coords) / N

    for d in layer:
        if d not in dots_for_remove:
            dots_for_remove.append(d.eat(layer, lsr * pararam, dots_for_remove))


    minus = 0
    for ddot in dots_for_remove:
        if ddot in layer:
            layer.remove(ddot)
            minus += 1


    layers.append(layer)



    return N - minus

def add_layer(N, prev_layer, lmax, pp_l):
    par = 0.9
    layer = []
    lsr = 0
    nd = None
    nd0 = None
    usd = N
    for i in range(0, N):
        usd = N
        angle = np.pi
        sch = len(layers)+1
        c00 = prev_layer[i].coords
        c0 = prev_layer[i].parentco
        '''
        while angle > np.pi * par and sch < 10:
            sch += 1
            v1 = prev_layer[(i - sch)%usd].coords - c0
            v2 = prev_layer[(i + sch) % usd].coords - c0
            angle = abs(np.arccos(np.dot(v1,v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
        '''
        vspdot = (prev_layer[(i - sch)%usd].coords + prev_layer[(i + sch) % usd].coords) / 2
        dir = (vspdot - c0) / np.linalg.norm(vspdot - c0)
        newdot = dir * 2 * min(lmax, 5 * (np.linalg.norm(prev_layer[(i - 1)%usd].coords - c0) + np.linalg.norm(prev_layer[(i + 1)%usd].coords - c0)) / 2)
        #newdot = dir * lmax
        newdot = c0 + newdot - proj(newdot, prev_layer[(i - sch)%usd].coords - prev_layer[(i + sch) % usd].coords)


        pd = nd
        nd = dot(newdot, prev_layer[i].coords)
        if i == 0:
            nd0 = nd
        layer.append(nd)
        if i != 0:
            triengles.append(triengle(pd, nd, prev_layer[(i - 1)%usd]))
            triengles.append(triengle(prev_layer[i], nd, prev_layer[(i - 1) % N]))






    triengles.append(triengle(nd, nd0, prev_layer[usd-1]))
    triengles.append(triengle(prev_layer[0], nd0, prev_layer[usd - 1]))

    lsr = 0
    dots_for_remove = []
    mid_dots = []
    prev_dot_wasnt_del = True
    pprev_dot_wasnt_del = True
    for d_num in range(len(layer)):
        lsr += np.linalg.norm(layer[d_num].coords - layer[(d_num + 1) % N].coords) / N

    for d in layer:
        if d not in dots_for_remove:
            dots_for_remove = dots_for_remove + d.eat(layer, lsr * pararam,dots_for_remove)


    '''
    for d_num in range(len(layer)):
        if np.linalg.norm(layer[d_num].coords - layer[(d_num + 1) % N].coords) < 0 * lsr and prev_dot_wasnt_del and pprev_dot_wasnt_del:
            print(1)
            pprev_dot_wasnt_del = prev_dot_wasnt_del
            prev_dot_wasnt_del = False
            mid_dot = dot((layer[d_num].coords + layer[(d_num + 1) % N].coords) / 2)
            mid_dots.append(mid_dot)
            d1 = layer[d_num]
            d2 = layer[(d_num + 1) % N]
            dots_for_remove.append(d1)
            dots_for_remove.append(d2)
            for t in d1.my_triengles + d2.my_triengles:
                if t in triengles:
                    triengles.remove(t)
            triengles.append(triengle(prev_layer[(d_num - 1) % N], prev_layer[(d_num) % N], mid_dot))
            triengles.append(triengle(prev_layer[(d_num - 1) % N], layer[(d_num - 1) % N], mid_dot))
            triengles.append(triengle(prev_layer[(d_num + 1) % N], prev_layer[(d_num) % N], mid_dot))
            triengles.append(triengle(prev_layer[(d_num + 1) % N], layer[(d_num + 2) % N], mid_dot))
        else:
            pprev_dot_wasnt_del = prev_dot_wasnt_del
            prev_dot_wasnt_del = True
        '''
    minus = 0
    for dddot in dots_for_remove:
        if dddot in layer:
            layer.remove(dddot)
            minus += 1


    layers.append(layer)

    print('N = ', N, '-', minus)
    return N - minus


def last_layer(N, previous_layer):
    center = np.array([0, 0, 0])
    for i in previous_layer:
        center = center + i.coords / N

    cd = dot(center)
    for i in range(len(previous_layer) - 1):

        triengles.append(triengle(previous_layer[i], previous_layer[i + 1], cd))
    triengles.append(triengle(previous_layer[0], previous_layer[N - 1], cd))








import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def draw_triangles(triangles):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    verts = []

    for triangle in triangles:
        v1 = triangle.d1.coords
        v2 = triangle.d2.coords
        v3 = triangle.d3.coords


        verts.append([v1, v2, v3])

    triangles_collection = Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='b', alpha=.25)

    t_values = np.linspace(0, 2 * np.pi, 100)
    #contour_points = contur(t_values)
    #ax.plot(contour_points[0], contour_points[1], contour_points[2], label='Contour', color='r')

    ax.add_collection3d(triangles_collection)
    #ax.scatter(verts)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([min(v[0] for tri in verts for v in tri) - 1, max(v[0] for tri in verts for v in tri) + 1])
    ax.set_ylim([min(v[1] for tri in verts for v in tri) - 1, max(v[1] for tri in verts for v in tri) + 1])
    ax.set_zlim([min(v[2] for tri in verts for v in tri) - 1, max(v[2] for tri in verts for v in tri) + 1])

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def draw_triangles2(triangles):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    verts = []
    X = []
    Y = []
    Z = []
    facecolors = []

    for triangle in triangles:
        v1 = triangle.d1.coords
        v2 = triangle.d2.coords
        v3 = triangle.d3.coords

        verts.append([v1, v2, v3])
        X.append(v1[0])
        X.append(v2[0])
        X.append(v3[0])
        Y.append(v1[1])
        Y.append(v2[1])
        Y.append(v3[1])
        Z.append(v1[2])
        Z.append(v2[2])
        Z.append(v3[2])



        # Используем среднюю Z-координату для вычисления цвета
        avg_z = (v1[2] + v2[2] + v3[2]) / 3
        normalized_color = 0.5  # Нормализация значения
        color = plt.cm.viridis(normalized_color)  # Используем colormap для получения цвета
        facecolors.append(color)

    min_z = min(v[2] for tri in verts for v in tri)
    max_z = max(v[2] for tri in verts for v in tri)

    triangles_collection = Poly3DCollection(verts, facecolors=facecolors, linewidths=1, edgecolors='b', alpha=.75)

    #ax.add_collection3d(triangles_collection)
    #ax.scatter(X, Y, Z, cmap='coolwarm', alpha=0.7)
    my_cmap = plt.get_cmap('hot')

    # Creating plot
    surf = ax.plot_trisurf(np.array(X), np.array(Y), np.array(Z),
                           cmap=my_cmap,
                           edgecolor='none')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([min(v[0] for tri in verts for v in tri) - 1, max(v[0] for tri in verts for v in tri) + 1])
    ax.set_ylim([min(v[1] for tri in verts for v in tri) - 1, max(v[1] for tri in verts for v in tri) + 1])
    ax.set_zlim([min_z - 1, max_z + 1])  # Используем min_z и max_z


# Пример использования функции
# triangles = [...]  # Определите ваши треугольники здесь
# draw_triangles(triangles)
plt.show()


def where_should_I_fall(layers, M):
    total_s = 0
    for lay in layers:
        for d in lay:
            char_l = max(np.dot(d.my_triengles[0].d1.coords - d.my_triengles[0].d2.coords, d.my_triengles[0].d1.coords - d.my_triengles[0].d2.coords)**0.5, np.dot(d.my_triengles[0].d1.coords - d.my_triengles[0].d3.coords, d.my_triengles[0].d1.coords - d.my_triengles[0].d3.coords)**0.5)
            char_l = char_l / M
            d.chrl = char_l
            sum_of_s0 = d.smt()
            total_s += sum_of_s0 / 3

            d.coords[0] += char_l
            dsx = d.smt()
            dsum_of_sx = abs(d.smt() - sum_of_s0)
            d.coords[0] -= char_l

            d.coords[1] += char_l
            dsy = d.smt()
            dsum_of_sy = abs(d.smt() - sum_of_s0)
            d.coords[1] -= char_l

            d.coords[2] += char_l
            dsz = d.smt()
            dsum_of_sz = abs(d.smt() - sum_of_s0)
            d.coords[2] -= char_l

            dir = np.array([dsum_of_sx, dsum_of_sy, dsum_of_sz])
            dir = dir / np.linalg.norm(dir)


            d.direction = -1 * dir
    return total_s

def moov_dots(layers):
    for lay in layers[1:]:
        for d in lay:
            if True:
                s0 = d.smt()
                d.coords = d.coords + d.direction * d.chrl
                s1 = d.smt()
                char_l = max(np.dot(d.my_triengles[0].d1.coords - d.my_triengles[0].d2.coords,
                                    d.my_triengles[0].d1.coords - d.my_triengles[0].d2.coords) ** 0.5,
                             np.dot(d.my_triengles[0].d1.coords - d.my_triengles[0].d3.coords,
                                    d.my_triengles[0].d1.coords - d.my_triengles[0].d3.coords) ** 0.5)
                char_l = char_l / M
                d.chrl = char_l
                if s1 > s0:
                    d.coords = d.coords - d.direction * 0.5 * d.chrl
                    s2 = d.smt()
                    if s2 > s0:
                        d.coords = d.coords - d.direction * 0.25 * d.chrl
                        s3 = d.smt()
                        if s3 > s0:
                            d.coords = d.coords - d.direction * 0.125 * d.chrl
                            s4 = d.smt()
                            if s4 > s0:
                                d.coords = d.coords - d.direction * 0.0625 * d.chrl
                                s5 = d.smt()
                                if s5 > s0:
                                    d.coords = d.coords - d.direction * 0.0625 * d.chrl
                                    d.moovable = False


adsad = add_layer1(n, layers[len(layers) - 1], LMAX, len(layers))

for _ in range(L - 1):
    print(adsad)
    adsad = add_layer(adsad, layers[len(layers) - 1], LMAX, layers[len(layers) - 2])





last_layer(adsad, layers[len(layers) - 1])
for t in triengles:
    if t.s() ==0:
        triengles.remove(t)

draw_triangles(triengles)


S0 = where_should_I_fall(layers, M)
for _ in range(MD):
    print(_, 'current S =', where_should_I_fall(layers, M))
    #where_should_I_fall(layers, M)
    moov_dots(layers)

S1 = where_should_I_fall(layers, M)

print(S0, '===>', S1)
draw_triangles2(triengles)
plt.show()













