import numpy as np

# rang de la cle p dans un dict
def get_rank(dict_p, p):
    for i, e in enumerate(dict_p):
        if e == p:
            return i
    return -1

# returns a list of coordinates from a shapely linestring, multilinestring or polygon 
def get_coords(shape):
    if shape.geom_type == 'MultiLineString':
        return shape[0].coords 
    elif shape.geom_type == 'Polygon':
        return shape.exterior.coords
    elif shape.geom_type == 'LineString':
        return shape.coords
    return []

class Squarer:
    """Initialize squaring object, with default weights and tolerance set in the constructor
    """
    def __init__(self, max_iter=1000, norm_tol=0.05,rtol=10, ftol=0.11, hrtol=10, pfixe=5, p90=100, p0=50, p45=10, switch_new=False):
        self.SWITCH_NEW = switch_new
        self.MAX_ITER = max_iter
        self.NORM_DIFF_TOL = norm_tol
        self.rightTol = rtol # 10 90° angles tolerance
        self.flatTol = ftol # 0.11 flat angles tolerance
        self.semiRightTol = hrtol # 7 45/135° angles tolerance
        self.poidsPtfFixe = pfixe #5
        self.poids90 = p90 #100
        self.poids0 = p0 #50
        self.poids45 = p45 #10

        self.point_shapes = {}
        self.lines_pindex = []
        self.indicesRight, self.indicesFlat, self.indicesHrAig, self.indicesHrObt = [], [], [], []

    # renvoie un dict dont la clé est le tuple d'un point, et la valeurs associee est
    # le tableau d'indices des lignes/poly dont le point fait partie
    def build_dict_of_unique_points(self, shapes):
        for i, s in enumerate(shapes):
            #coords = s[0].coords if shapes[0].geom_type == 'MultiLineString' else s.exterior.coords
            coords = get_coords(s)
            for p in coords:
                #print(p) # tuple
                if p in self.point_shapes:
                    if i not in self.point_shapes[p]:
                        self.point_shapes[p].append(i)
                else:
                    self.point_shapes[p] = [i]


    # renvoie une liste où chaque ligne est constituée d'une liste des indices des points
    # composant la ligne/poly à cet indice 
    def build_pindex_for_shapes(self, shapes):
        for s in shapes:
            index_points = []
            #coords = s[0].coords if shapes[0].geom_type == 'MultiLineString' else s.exterior.coords
            coords = get_coords(s)
            for p in coords:
                index_points.append(get_rank(self.point_shapes, p))
            self.lines_pindex.append(index_points)

    # renvoie le rang du point d'indice idx_p dans la shape d'indice idx_s
    def get_rank_point_in_shape(self, idx_p, idx_s):
        for i, idx_pp in enumerate(self.lines_pindex[idx_s]):
            if idx_pp == idx_p: return i
        return -1

    # renvoie la liste des angles potentiels autour du point d'indice idx_p sous forme de 
    # triplets d'indices [[idx_prec, idx_p, idx_suiv]...]
    def get_angle_triplets(self, idx_p, unik_points):
        p = unik_points[idx_p]
        lines_containing_p = self.point_shapes[p]
        if len(lines_containing_p) == 1:
            idx_l = lines_containing_p[0]
            r = self.get_rank_point_in_shape(idx_p, idx_l)
            # single node, no angle
            if r == 0 or r == len(self.lines_pindex[idx_l]) - 1 : 
                return []
            # interior point of a line
            idx_prec = self.lines_pindex[idx_l][r - 1]
            idx_suiv = self.lines_pindex[idx_l][r + 1]
            #print(f'POINT({p[0]} {p[1]})')
            return [[idx_prec, idx_p, idx_suiv]]
        # points intersecting multiple lines
        else:
            triplets = []
            #print(f'POINT({p[0]} {p[1]})')
            for i in range(len(lines_containing_p) - 1):
                idx_l1 = lines_containing_p[i]
                r = self.get_rank_point_in_shape(idx_p, idx_l1)
                idx_prec = self.lines_pindex[idx_l1][r - 1] if r > 0 else self.lines_pindex[idx_l1][r + 1]
                for j in range(i + 1, len(lines_containing_p)):
                    idx_l2 = lines_containing_p[j]
                    r = self.get_rank_point_in_shape(idx_p, idx_l2)
                    idx_suiv = self.lines_pindex[idx_l2][r - 1] if r > 0 else self.lines_pindex[idx_l2][r + 1]
                    triplets.append([idx_prec, idx_p, idx_suiv])
            return triplets

    def get_vecs_around(self, t, unik_points) : # t = [idx_prec, idx_p, idx_suiv]
        """ return vectors formed by a triplet of indexes
        """
        pr, p, s = unik_points[t[0]], unik_points[t[1]], unik_points[t[2]]
        if self.SWITCH_NEW:
            v1 = np.array([pr[0] - p[0], pr[1] - p[1]])
        else:
            v1 = np.array([p[0] - pr[0], p[1] - pr[1]])
        v2 = np.array([s[0] - p[0], s[1] - p[1]])
        return v1, v2

    def idx_angles_remarquables(self, unik_points):
        #unik_points = list(self.point_shapes)
        rTol = np.cos((np.pi / 2) - self.rightTol * np.pi / 180)
        hrTol1 = np.cos((np.pi / 4) - self.semiRightTol * np.pi / 180)
        hrTol2 = np.cos((np.pi / 4) + self.semiRightTol * np.pi / 180)
        for idx_p in range(len(unik_points)):
            triplets = self.get_angle_triplets(idx_p, unik_points)
            for t in triplets:
                v1, v2 = self.get_vecs_around(t, unik_points)
                n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
                v1n = v1 / n1 if n1 != 0. else np.array([0.,0.]) #n1
                v2n = v2 / n2 if n2 != 0. else np.array([0.,0.]) #n2
                dot = v1n.dot(v2n)
                cross = np.cross(v1n, v2n).item(0)
                if (np.abs(dot) <= rTol):
                    self.indicesRight.append(t)
                elif (cross <= self.flatTol):
                    self.indicesFlat.append(t)
                elif (dot <= hrTol1 and dot >= hrTol2):
                    self.indicesHrAig.append(t)
                elif (dot >= -hrTol1 and dot <= -hrTol2):
                    self.indicesHrObt.append(t)
        print(f'potential angles -- R: {len(self.indicesRight)} - F: {len(self.indicesFlat)} - HRa: {len(self.indicesHrAig)} - HRo: {len(self.indicesHrObt)}')


    def get_Y(self, unik_points):
        """ Observation vector
        """
        nb_points = len(unik_points)
        self.Y = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrObt) + len(self.indicesHrAig))
        for i, p in enumerate(unik_points):
            self.Y[2*i] = p[0]
            self.Y[2*i+1] = p[1]
        offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat)
        for i, t in enumerate(self.indicesHrAig):
            v1, v2 = self.get_vecs_around(t, unik_points)
            d = np.linalg.norm(v1) * np.linalg.norm(v2) * np.cos(np.pi / 4)
            self.Y[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrAig)
        for i, t in enumerate(self.indicesHrObt):
            v1, v2 = self.get_vecs_around(t, unik_points)
            d = np.linalg.norm(v1) * np.linalg.norm(v2) * np.cos(3 * np.pi / 4)
            self.Y[offset + i] = d

    # B = Y - S(Xcourant)
    def get_B(self, points):
        nb_points = len(points)
        S = np.zeros(2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrObt) + len(self.indicesHrAig))
        for i, p in enumerate(points):
            S[2*i] = p[0]
            S[2*i+1] = p[1]
        offset = 2 * nb_points
        for i, t in enumerate(self.indicesRight):
            v1, v2 = self.get_vecs_around(t, points)
            d = v1.dot(v2) 
            S[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight)
        for i, t in enumerate(self.indicesFlat):
            v1, v2 = self.get_vecs_around(t, points)
            d = np.cross(v1, v2).item(0) 
            S[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat)
        for i, t in enumerate(self.indicesHrAig):
            v1, v2 = self.get_vecs_around(t, points)
            d = v1.dot(v2) 
            S[offset + i] = d
        offset = 2 * nb_points + len(self.indicesRight) + len(self.indicesFlat) + len(self.indicesHrAig)
        for i, t in enumerate(self.indicesHrObt):
            v1, v2 = self.get_vecs_around(t, points)
            d = v1.dot(v2) 
            S[offset + i] = d
        return S - self.Y

    # Weight Matrix
    # n = 2 * nb_points + indicesRight.size() + indicesFlat.size() + indicesHrAig.size() + indicesHrObt.size()
    def get_P(self):
        nb_points = len(self.point_shapes)
        nb_rights, nb_flats =  len(self.indicesRight), len(self.indicesFlat)
        nb_half_rights = len(self.indicesHrAig) + len(self.indicesHrObt)
        wfix = np.full(2*nb_points, self.poidsPtfFixe)
        wRight = np.full(nb_rights, self.poids90)
        wFlat = np.full(nb_flats, self.poids0)
        wHr = np.full(nb_half_rights, self.poids45)
        self.P = np.diag(np.concatenate((wfix, wRight, wFlat, wHr)))

    # old way
    def partial_derivatives_dotpo(self, points, indices):
        nb_points = len(points)
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = p[0] - s[0]
            dfy = p[1] - s[1]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = s[0] - 2*p[0] + pr[0]
            dfy = s[1] - 2*p[1] + pr[1]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = p[0] - pr[0]
            dfy = p[1] - pr[1]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m

    def partial_derivatives_crosso(self, points, indices):
        nb_points = len(points) #- 1
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = p[1] - s[1]
            dfy = -p[0] + s[0]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = s[1] - pr[1]
            dfy = -s[0] + pr[0]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = -p[1] + pr[1]
            dfy = p[0] - pr[0]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m

    ## new vectors
    def partial_derivatives_dotpn(self, points, indices):
        nb_points = len(points)
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = s[0] - p[0]
            dfy = s[1] - p[1]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = -pr[0] + 2*p[0] - s[0]
            dfy = -pr[1] + 2*p[1] - pr[1]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = pr[0] - p[0]
            dfy = pr[1] - p[1]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m

    def partial_derivatives_crossn(self, points, indices):
        nb_points = len(points) #- 1
        nb_indices = len(indices)
        m = np.zeros((nb_indices, 2*nb_points))
        for i, t in enumerate(indices):
            idx_prec, idx, idx_suiv = t[0], t[1], t[2]
            pr, p, s = points[t[0]], points[t[1]], points[t[2]]
            # df en Xi-1, Yi-1
            dfx = s[1] - p[1]
            dfy = -s[0] + p[0]
            m[i][2*idx_prec] = dfx
            m[i][2*idx_prec + 1] = dfy
            # df en Xi, Yi
            dfx = -s[1] + pr[1]
            dfy = s[0] - pr[0]
            m[i][2*idx] = dfx
            m[i][2*idx + 1] = dfy
            # df en Xi+1, Yi+1
            dfx = -pr[1] + p[1]
            dfy = pr[0] - p[0]
            m[i][2*idx_suiv] = dfx
            m[i][2*idx_suiv + 1] = dfy
        return m

    def partial_derivatives_dotp(self, points, indices): 
        if self.SWITCH_NEW:
            return self.partial_derivatives_dotpn(points, indices) 
        else: 
            return self.partial_derivatives_dotpo(points, indices)
    
    def partial_derivatives_cross(self, points, indices): 
        if self.SWITCH_NEW:
            return self.partial_derivatives_crossn(points, indices) 
        else: 
            return self.partial_derivatives_crosso(points, indices)

    def get_A(self, points):
        nb_points = len(points) #- 1
        id = np.identity(2 * nb_points)
        partialR = self.partial_derivatives_dotp(points, self.indicesRight)
        partialCross = self.partial_derivatives_cross(points, self.indicesFlat)
        partialHr1 = self.partial_derivatives_dotp(points, self.indicesHrAig)
        partialHr2 = self.partial_derivatives_dotp(points, self.indicesHrObt)
        a = np.vstack((id, partialR))
        a = np.vstack((a, partialCross))
        a = np.vstack((a, partialHr1))
        a = np.vstack((a, partialHr2))
        return a

    def compute_dx(self, points):
        A = self.get_A(points)
        B = self.get_B(points)
        atp = A.T @ self.P 
        atpa = atp @ A
        atpb = atp @ B
        #dx = np.linalg.lstsq(atpa, atpb)
        #dx = np.linalg.inv(atpa) @ atpb
        dx = np.linalg.solve(atpa, atpb)
        return dx
    
    def prepare_square(self, shapes):
        if len(shapes) == 0:
            return np.array([])
        self.geom_type = shapes[0].geom_type
        self.build_dict_of_unique_points(shapes)
        self.build_pindex_for_shapes(shapes)
        unik_points = list(self.point_shapes)
        self.idx_angles_remarquables(unik_points)
        self.get_Y(unik_points)
        self.get_P()
        return np.array(unik_points)
    
    def square(self, shapes):
        """squares a collection of shapely multilinestrings or polygons
        returns a numpy array of the points after the least square process
        """
        points = self.prepare_square(shapes)
        nb_points = len(points)
        for i in range(self.MAX_ITER):
            dx = self.compute_dx(points)
            points -= dx.reshape((nb_points, 2))
            print(i, np.linalg.norm(dx, ord=np.inf))
            if np.linalg.norm(dx, ord=np.inf) < self.NORM_DIFF_TOL:
                break
        self.nb_iters = i
        return points


    # rebuild shapes with updated points from least square process
    def get_shapes_from_new_points(self, original_shapes, new_points):
        """rebuild a list of coordinates from the original collection  of shapes
        and the points obtained from the square process
        """
        unik_points = list(self.point_shapes)
        new_s = []
        for l in original_shapes:
            #coords, is_poly = (l[0].coords, False) if l.geom_type == 'MultiLineString' else (l.exterior.coords, True)
            coords = get_coords(l)
            is_poly = True if self.geom_type == 'Polygon' else False
            size = len(coords)
            new_s.append(np.zeros((size, 2)) )
        for idx_p, p in enumerate(unik_points):
            index_of_lines = self.point_shapes[p] # point_lines_idx[p]
            #print(index_of_lines)
            for idx_l in index_of_lines:
                r = self.get_rank_point_in_shape(idx_p, idx_l)
                #print(idx_l, r, new_points[idx_p])
                new_s[idx_l][r] = np.array(new_points[idx_p])
                if r == 0 and is_poly:
                    new_s[idx_l][-1] = np.array(new_points[idx_p])
        return new_s



if __name__ == '__main__':
    import utils
    from shapely.wkt import loads
    from shapely.geometry import Polygon, LineString

    linesf = 'clermont_bbox.geojson'
    linesf = 'clermont_small.geojson'
    lines = utils.load_geojson(linesf)
    wkt = "POLYGON((342763.9 7685274.5, 342763.9 7685274.1, 342767.0 7685273.9, 342767.1 7685274.4, 342770.9 7685274.1, 342770.8 7685272.4, 342772.1 7685272.3, 342772.0 7685269.0, 342774.9 7685268.9, 342775.1 7685272.1, 342776.3 7685272.1, 342776.3 7685273.9, 342779.6 7685273.9, 342779.7 7685273.5, 342782.7 7685273.4, 342782.8 7685273.8, 342785.9 7685273.7, 342785.9 7685271.9, 342787.3 7685271.8, 342787.2 7685268.4, 342790.2 7685268.2, 342790.3 7685271.7, 342791.4 7685271.7, 342791.5 7685273.4, 342794.7 7685273.4, 342794.7 7685272.8 , 342797.7 7685272.7 , 342797.8 7685273.1 , 342801.5 7685273.1 , 342801.7 7685271.4 , 342803.1 7685271.5 , 342803.6 7685268.2 , 342806.4 7685268.6 , 342805.9 7685272.0 , 342807.3 7685272.3 , 342806.9 7685273.9 , 342810.6 7685275.1 , 342810.8 7685274.8 , 342813.6 7685275.7 , 342813.5 7685276.1 , 342817.2 7685277.1 , 342817.9 7685275.8 , 342819.1 7685276.5 , 342820.7 7685273.6 , 342823.3 7685275.1 , 342821.6 7685278.0 , 342822.7 7685278.7 , 342821.7 7685280.1 , 342824.6 7685282.4 , 342825.0 7685282.0 , 342827.3 7685284.0 , 342827.0 7685284.5 , 342829.8 7685286.9 , 342831.2 7685285.7 , 342832.0 7685286.5 , 342834.6 7685284.5 , 342841.3 7685292.7 , 342821.1 7685309.4 , 342819.1 7685306.7 , 342817.6 7685307.8 , 342815.7 7685309.1 , 342815.1 7685308.3 , 342813.2 7685305.8 , 342814.6 7685304.5 , 342816.0 7685303.2 , 342818.3 7685301.4 , 342817.0 7685299.7 , 342821.7 7685295.3 , 342823.7 7685293.5 , 342823.5 7685293.3 , 342823.4 7685293.2 , 342823.3 7685293.1 , 342823.1 7685293.0 , 342822.9 7685292.8 , 342822.8 7685292.7 , 342822.7 7685292.6 , 342822.6 7685292.5 , 342822.5 7685292.4 , 342822.4 7685292.2 , 342822.2 7685292.1 , 342822.0 7685291.9 , 342821.9 7685291.8 , 342821.8 7685291.7 , 342821.6 7685291.5 , 342821.4 7685291.3 , 342821.3 7685291.2 , 342819.3 7685293.2 , 342817.2 7685291.6 , 342818.8 7685289.4 , 342818.5 7685289.2 , 342818.3 7685289.1 , 342818.1 7685288.9 , 342817.9 7685288.7 , 342817.7 7685288.6 , 342817.5 7685288.4 , 342817.3 7685288.3 , 342817.1 7685288.1 , 342816.9 7685288.0 , 342816.7 7685287.8 , 342816.5 7685287.7 , 342816.2 7685287.5 , 342814.8 7685289.8 , 342812.4 7685288.3 , 342813.5 7685286.0 , 342813.2 7685285.8 , 342813.0 7685285.7 , 342812.7 7685285.6 , 342812.5 7685285.5 , 342812.2 7685285.3 , 342811.8 7685285.1 , 342811.5 7685285.0 , 342811.0 7685284.8 , 342810.8 7685284.7 , 342810.5 7685284.5 , 342809.4 7685287.1 , 342806.6 7685286.2 , 342807.3 7685283.6 , 342807.1 7685283.6 , 342806.8 7685283.5 , 342806.5 7685283.5 , 342806.2 7685283.4 , 342806.0 7685283.3 , 342805.7 7685283.3 , 342805.4 7685283.2 , 342805.2 7685283.1 , 342804.9 7685283.1 , 342804.6 7685283.0 , 342804.4 7685283.0 , 342804.1 7685282.9 , 342803.6 7685285.6 , 342800.9 7685285.1 , 342801.2 7685282.6 , 342801.0 7685282.6 , 342800.7 7685282.6 , 342800.4 7685282.6 , 342800.2 7685282.5 , 342799.9 7685282.5 , 342799.7 7685282.5 , 342799.4 7685282.5 , 342799.1 7685282.4 , 342798.9 7685282.4 , 342798.6 7685282.4 , 342798.4 7685282.4 , 342798.1 7685282.3 , 342797.8 7685282.3 , 342797.8 7685284.8 , 342795.1 7685284.6 , 342794.9 7685282.4 , 342794.5 7685282.4 , 342794.3 7685282.4 , 342794.0 7685282.4 , 342793.7 7685282.4 , 342793.4 7685282.4 , 342793.1 7685282.4 , 342792.5 7685282.4 , 342791.9 7685282.4 , 342791.3 7685282.4 , 342791.0 7685282.4 , 342790.7 7685282.5 , 342790.5 7685282.5 , 342790.6 7685284.7 , 342790.6 7685285.0 , 342790.2 7685285.0 , 342789.8 7685285.1 , 342789.3 7685285.1 , 342788.9 7685285.1 , 342788.5 7685285.2 , 342788.0 7685285.2 , 342787.7 7685285.2 , 342787.7 7685284.7 , 342787.6 7685282.6 , 342787.4 7685282.6 , 342787.1 7685282.7 , 342786.8 7685282.7 , 342786.5 7685282.7 , 342786.2 7685282.7 , 342785.8 7685282.7 , 342785.5 7685282.8 , 342785.2 7685282.8 , 342784.9 7685282.8 , 342784.6 7685282.8 , 342783.0 7685282.8 , 342783.0 7685284.9 , 342783.0 7685285.4 , 342782.8 7685285.4 , 342782.6 7685285.4 , 342782.2 7685285.4 , 342781.9 7685285.4 , 342781.7 7685285.4 , 342781.5 7685285.5 , 342781.4 7685285.5 , 342781.2 7685285.5 , 342781.0 7685285.5 , 342780.8 7685285.5 , 342780.7 7685285.5 , 342780.5 7685285.5 , 342780.3 7685285.6 , 342780.0 7685285.6 , 342780.0 7685285.1 , 342779.9 7685282.9 , 342775.4 7685283.1 , 342775.4 7685285.2 , 342775.4 7685285.7 , 342772.4 7685285.8 , 342772.4 7685285.3 , 342772.4 7685283.2 , 342767.7 7685283.4 , 342767.7 7685285.6 , 342767.7 7685286.0 , 342764.9 7685286.2 , 342764.5 7685283.6 , 342761.5 7685284.0 , 342762.1 7685286.5 , 342759.2 7685287.3 , 342758.3 7685284.8 , 342758.1 7685284.9 , 342757.8 7685285.0 , 342757.5 7685285.1 , 342757.1 7685285.2 , 342756.8 7685285.4 , 342756.5 7685285.5 , 342756.2 7685285.7 , 342755.8 7685285.8 , 342755.5 7685286.0 , 342756.6 7685288.3 , 342756.4 7685288.4 , 342756.2 7685288.5 , 342756.1 7685288.6 , 342756.0 7685288.6 , 342755.8 7685288.7 , 342755.6 7685288.8 , 342755.4 7685288.9 , 342755.3 7685289.0 , 342755.1 7685289.1 , 342754.8 7685289.2 , 342754.7 7685289.3 , 342754.6 7685289.3 , 342754.4 7685289.4 , 342754.2 7685289.5 , 342754.0 7685289.6 , 342752.7 7685287.4 , 342752.4 7685287.6 , 342752.2 7685287.8 , 342751.9 7685287.9 , 342751.6 7685288.1 , 342751.3 7685288.3 , 342750.8 7685288.6 , 342750.6 7685288.8 , 342750.3 7685289.0 , 342750.0 7685289.2 , 342751.5 7685291.3 , 342749.2 7685293.1 , 342748.9 7685292.7 , 342747.5 7685291.2 , 342747.3 7685291.4 , 342747.0 7685291.6 , 342746.8 7685291.8 , 342746.3 7685292.2 , 342746.0 7685292.4 , 342745.8 7685292.6 , 342745.5 7685292.8 , 342745.3 7685293.0 , 342745.0 7685293.2 , 342744.7 7685293.4 , 342744.5 7685293.6 , 342744.2 7685293.8 , 342744.1 7685294.0 , 342745.4 7685295.6 , 342745.8 7685296.0 , 342743.5 7685297.9 , 342741.3 7685299.8 , 342741.0 7685299.5 , 342739.6 7685297.7 , 342736.1 7685300.6 , 342737.5 7685302.3 , 342737.8 7685302.7 , 342735.5 7685304.6 , 342735.2 7685304.3 , 342733.8 7685302.5 , 342730.2 7685305.5 , 342731.7 7685307.2 , 342731.9 7685307.5 , 342729.7 7685309.4 , 342729.4 7685309.0 , 342728.0 7685307.4 , 342724.5 7685310.3 , 342725.9 7685311.9 , 342726.2 7685312.3 , 342723.9 7685314.3 , 342723.5 7685313.8 , 342722.2 7685312.2 , 342718.6 7685315.2 , 342712.6 7685308.0 , 342715.2 7685305.8 , 342714.2 7685304.7 , 342715.3 7685303.8 , 342713.1 7685301.2 , 342715.4 7685299.3 , 342717.6 7685301.9 , 342718.5 7685301.1 , 342719.6 7685302.5 , 342722.3 7685300.4 , 342721.9 7685300.0 , 342724.3 7685298.0 , 342724.5 7685298.3 , 342727.0 7685296.3 , 342725.9 7685294.9 , 342726.8 7685294.1 , 342724.9 7685291.7 , 342727.2 7685289.8 , 342727.6 7685290.2 , 342729.2 7685292.1 , 342730.0 7685291.3 , 342731.2 7685292.7 , 342733.8 7685290.7 , 342733.5 7685290.3 , 342738.0 7685286.6 , 342738.2 7685287.0 , 342740.8 7685284.9 , 342739.5 7685283.1 , 342740.3 7685282.4 , 342738.6 7685280.4 , 342738.5 7685280.2 , 342740.7 7685278.3 , 342742.8 7685280.9 , 342744.0 7685279.9 , 342745.0 7685281.3 , 342748.7 7685279.5 , 342748.5 7685279.0 , 342751.1 7685277.7 , 342751.3 7685278.0 , 342755.1 7685276.0 , 342754.7 7685274.3 , 342755.9 7685273.9 , 342755.1 7685271.0 , 342757.9 7685270.2 , 342758.8 7685273.2 , 342760.1 7685273.0 , 342760.3 7685274.8 , 342763.9 7685274.5))"
    #lines = [loads(wkt)]
    sq = Squarer(pfixe=5)
    with np.printoptions(precision=3, suppress=True):
        points = sq.square(lines)

    new_s = sq.get_shapes_from_new_points(lines, points)
    is_line = True if lines[0].geom_type == 'MultiLineString' else False
    for e in new_s:
        if is_line:
            print(LineString(e))
        else:
            print(Polygon(e)) # LineString(e))

    print("done in", sq.nb_iters, "iters")