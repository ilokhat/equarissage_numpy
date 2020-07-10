from lsq.squarer import Squarer
from lsq.utils import load_geojson, build_wkts_from_coords, load_shapefile

geojson = 'data/clermont_small.geojson'
lines = load_geojson(geojson)

sq = Squarer(pfixe=5000)
new_points = sq.square(lines)
new_coords = sq.get_shapes_from_new_points(lines, new_points)
res = build_wkts_from_coords(new_coords)
print("done in", sq.nb_iters, "iters")
for l in res:
    print(l)

# polys = load_shapefile(shapefile)
# shapefile = 'buildings.shp'
# res = []
# for p in polys:
#     sq = Squarer(pfixe=5,switch_new=True)
#     new_points = sq.square([p])
#     new_coords = sq.get_shapes_from_new_points([p], new_points)
#     res.append(build_wkts_from_coords(new_coords, shape_type='Polygon')[0])

# for p in res:
#     print(p)
