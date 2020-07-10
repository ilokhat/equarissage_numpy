from shapely.geometry import shape, GeometryCollection, Polygon, LineString
import json
import fiona

def load_geojson(geojson):
    """Load a geojson file as a shapely geometry collection
    (expects a geojson file containing multilinestrings or polygons)
    """
    with open(geojson) as f:
        features = json.load(f)["features"]
    coll = GeometryCollection([shape(feature["geometry"]) for feature in features])   
    return coll

def load_shapefile(shapefile):
    """Load a geojson file as a geometry collection
    (expects a multilinestring or polygon collection)
    """
    with fiona.open(shapefile) as shp:
        feats = []
        for f in shp:
            feats.append(shape(f['geometry']))
        return GeometryCollection(feats)

def build_wkts_from_coords(coords_list, shape_type='MultiLineString'):
    """returns a list of wkts from from a list of coordinates and their type
    ('MultiLineString' by default, or 'Polygon')
    """
    shapes = []
    for e in coords_list:
        if shape_type == 'MultiLineString':
            shapes.append(LineString(e))
        else:
            shapes.append(Polygon(e))
    return shapes