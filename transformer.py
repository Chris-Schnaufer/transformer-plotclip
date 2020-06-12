"""Performs plot-level clipping on geo referenced files
"""

import argparse
import copy
import datetime
import json
import logging
import os
import re
import subprocess
from typing import Optional
import requests
from osgeo import gdal, ogr, osr
import yaml
import liblas
import numpy as np
from dbfread import DBF

import configuration
import transformer_class

BETYDB_URL = "https://terraref.ncsa.illinois.edu/bety"
BETYDB_LOCAL_CACHE_FOLDER = os.environ.get('BETYDB_LOCAL_CACHE_FOLDER', '/home/extractor/')

BETYDB_CULTIVARS = None
BETYDB_TRAITS = None
BETYDB_EXPERIMENTS = None


class __internal__():
    """Class for internal use only functions
    """

    def __init__(self):
        """Initializes class instance
        """

    @staticmethod
    def get_bety_key():
        """return key from environment or ~/.betykey if it exists."""

        key = os.environ.get('BETYDB_KEY', '')
        if key:
            return key

        keyfile_path = os.path.expanduser('~/.betykey')
        if os.path.exists(keyfile_path):
            keyfile = open(keyfile_path, "r")
            return keyfile.readline().strip()

        raise RuntimeError("BETYDB_KEY not found. Set environmental variable " +
                           "or create $HOME/.betykey.")

    @staticmethod
    def get_bety_url(path=''):
        """return betydb url from environment with optional path

        Of 3 options string join, os.path.join and urlparse.urljoin, os.path.join
        is the best at handling excessive / characters.
        """

        url = os.environ.get('BETYDB_URL', BETYDB_URL)
        return os.path.join(url, path)

    @staticmethod
    def get_bety_api(endpoint=None):
        """return betydb API based on betydb url"""

        url = __internal__.get_bety_url(path='api/v1/{}'.format(endpoint))
        return url

    @staticmethod
    def query(endpoint="search", **kwargs):
        """return betydb API results.

        This is general function for querying the betyDB API. It automatically
        decodes the json response if one is returned.
        """

        payload = {'key': __internal__.get_bety_key()}
        payload.update(kwargs)

        req = requests.get(__internal__.get_bety_api(endpoint), params=payload)
        req.raise_for_status()
        return req.json()

    @staticmethod
    def get_experiments(**kwargs):
        """Return cleaned up array from query() for the experiments table.
            If global variable isn't populated, check if a local file is present and read from it if so.
            This is for deployments where data is pre-fetched (e.g. for a Condor job).
            Otherwise the BETY API will be called.
            In either case, data will be kept in memory for subsequent calls.
        """
        global BETYDB_EXPERIMENTS

        if BETYDB_EXPERIMENTS is None:
            cache_file = os.path.join(BETYDB_LOCAL_CACHE_FOLDER, "bety_experiments.json")
            if os.path.exists(cache_file):
                with open(cache_file) as infile:
                    query_data = json.load(infile)
                    if query_data:
                        if 'associations_mode' in kwargs:
                            BETYDB_EXPERIMENTS = query_data
                        return [t["experiment"] for t in query_data['data']]
            else:
                query_data = __internal__.query(endpoint="experiments", **kwargs)
                if query_data:
                    if 'associations_mode' in kwargs:
                        BETYDB_EXPERIMENTS = query_data
                    return [t["experiment"] for t in query_data['data']]
        else:
            return [t["experiment"] for t in BETYDB_EXPERIMENTS['data']]

        return []

    @staticmethod
    def get_sites(filter_date='', include_halves=False, **kwargs):
        """Return a site array from query() from the sites table.

        e.g.
                get_sites(city="Maricopa")
                get_sites(sitename="MAC Field Scanner Season 4 Range 4 Column 6")
                get_sites(contains="-111.97496613200647,33.074671230742446")

          filter_date -- YYYY-MM-DD to filter sites to specific experiment by date
        """

        if not filter_date:
            # SCENARIO I - NO FILTER DATE
            # Basic query, efficient even with 'containing' parameter.
            query_data = __internal__.query(endpoint="sites", limit='none', **kwargs)
            if query_data:
                return [t["site"] for t in query_data['data']]
        else:
            # SCENARIO II - YES FILTER DATE
            # Get experiments by date and return all associated sites, optionally filtering by location.
            targ_date = datetime.datetime.strptime(filter_date, '%Y-%m-%d')
            query_data = __internal__.get_experiments(associations_mode='full_info', limit='none', **kwargs)
            if query_data:
                results = []
                for exp in query_data:
                    start = datetime.datetime.strptime(exp['start_date'], '%Y-%m-%d')
                    end = datetime.datetime.strptime(exp['end_date'], '%Y-%m-%d')
                    if start <= targ_date <= end and 'sites' in exp:
                        for one_entry in exp['sites']:
                            site = one_entry['site']
                            if (site["sitename"].endswith(" W") or site["sitename"].endswith(" E")) \
                                    and not include_halves:
                                continue
                            if 'containing' in kwargs:
                                # Need to filter additionally by geometry
                                site_geom = ogr.CreateGeometryFromWkt(site['geometry'])
                                coords = kwargs['containing'].split(",")
                                pt_geom = ogr.CreateGeometryFromWkt("POINT(%s %s)" % (coords[1], coords[0]))
                                if site_geom.Intersects(pt_geom):
                                    if site not in results:
                                        results.append(site)
                            else:
                                # If no containing parameter, include all sites
                                if site not in results:
                                    results.append(site)
                return results
            logging.error("No experiment data could be retrieved.")

        return []

    @staticmethod
    def get_site_boundaries(filter_date='', **kwargs):
        """Get a dictionary of site GeoJSON bounding boxes filtered by standard arguments.

        filter_date -- YYYY-MM-DD to filter sites to specific experiment by date

        Returns:
            {
                'sitename_1': 'geojson bbox',
                'sitename_2': 'geojson bbox',
                ...
             }
        """

        sitelist = __internal__.get_sites(filter_date, **kwargs)
        bboxes = {}

        for site in sitelist:
            geom = ogr.CreateGeometryFromWkt(site['geometry'])

            if geom:
                bboxes[site['sitename']] = __internal__.geometry_to_geojson(geom, 'EPSG', '4326')
            else:
                logging.error("Site boundary geometry is invalid for site: %s", site['sitename'])

        return bboxes

    @staticmethod
    def load_shapefile(shapefile: str, plot_column_name: str = None) -> list:
        """Loads the shapefile and returns a list of geometries
        Arguments:
            shapefile: the path of the shapefile to load
            plot_column_name: the name of the column containing the plot names
        Return:
            Returns the list of loaded plot boundaries
        """
        plots = []

        # Setup to parse the shapefile
        shape_in = ogr.Open(shapefile)
        layer_name = os.path.split(os.path.splitext(shapefile)[0])[1]
        if isinstance(layer_name, (bytes, bytearray)):
            layer_name = layer_name.decode('utf8')
        layer = shape_in.GetLayer(layer_name)
        feature = layer.GetNextFeature()
        layer_sr = layer.GetSpatialRef()

        # Get out target Spatial Reference
        default_sr = osr.SpatialReference()
        default_sr.ImportFromEPSG(4326)

        # Get the column name to use, if possible
        dbffile = os.path.splitext(shapefile)[0] + ".dbf"
        shape_table = DBF(dbffile, lowernames=True, ignore_missing_memofile=True)
        shape_rows = iter(list(shape_table))

        # Make sure if we have the column name of plot-names specified that it exists in
        # the shapefile
        column_names = shape_table.field_names
        if plot_column_name is not None:
            if plot_column_name not in column_names:
                raise RuntimeError("Shapefile does not have specified plot name column '%s'" % plot_column_name)

        # Lookup a plot name field to use
        if plot_column_name is None:
            for one_name in column_names:
                if one_name == "observationUnitName":
                    plot_column_name = one_name
                    break
                elif (one_name.find('plot') >= 0) and ((one_name.find('name') >= 0) or one_name.find('id')):
                    plot_column_name = one_name
                    break
                elif one_name == 'id':
                    plot_column_name = one_name
                    break
        if plot_column_name is None:
            logging.warning("Shapefile data does not have a known plot name field; using default naming")

        # Loop through each polygon and extract plot level data
        alternate_plot_id = 0
        logging.warning("About to loop through features")
        while feature:
            # Current geometry to extract
            plot_poly = feature.GetGeometryRef()
            if layer_sr:
                plot_poly.AssignSpatialReference(layer_sr)
            plot_spatial_ref = plot_poly.GetSpatialReference()

            # Check if we need to convert the polygon coordinate system
            if not plot_spatial_ref.IsSame(default_sr):
                plot_poly = __internal__.convert_geometry(plot_poly, default_sr)

            # Determine the plot name to use
            plot_name = None
            alternate_plot_id = alternate_plot_id + 1
            if shape_rows and plot_column_name:
                try:
                    row = next(shape_rows)
                    if plot_column_name in row:
                        plot_name = str(row[plot_column_name])
                except StopIteration:
                    pass
            if not plot_name:
                plot_name = "plot_" + str(alternate_plot_id)

            plots[plot_name] = __internal__.geometry_to_geojson(plot_poly)

        return plots

    @staticmethod
    def convert_json_geometry(geojson: str, new_spatialreference: Optional[osr.SpatialReference]) -> str:
        """Converts geojson geometry to new spatial reference system and
           returns the new geometry as geojson.
        Arguments:
            geojson: The geometry to transform
            new_spatialreference: The spatial reference to change to
        Return:
            The transformed geometry as geojson or the original geojson. If
            either the new Spatial Reference parameter is None, or the geojson
            doesn't have a spatial reference, or the geojson isn't a
            valid geometry, then the original geojson is returned.
        """
        if not new_spatialreference or not geojson:
            return geojson

        geom_yaml = yaml.safe_load(geojson)
        geometry = ogr.CreateGeometryFromJson(json.dumps(geom_yaml))

        if not geometry:
            return geojson

        new_geometry = __internal__.convert_geometry(geometry, new_spatialreference)
        if not new_geometry:
            return geojson

        try:
            return __internal__.geometry_to_geojson(new_geometry)
        except Exception as ex:
            logging.warning("Exception caught while transforming geojson: %s", str(ex))
            logging.warning("    Returning original geojson")

        return geojson

    @staticmethod
    def geometry_to_geojson(geom: ogr.Geometry, alt_coord_type: str = None, alt_coord_code: str = None) -> str:
        """Converts a geometry to geojson.
        Args:
            geom: The geometry to convert to JSON
            alt_coord_type: the alternate geographic coordinate system type if geometry doesn't have one defined
            alt_coord_code: the alternate geographic coordinate system associated with the type
        Returns:
            The geojson string for the geometry
        Note:
            If the geometry doesn't have a spatial reference associated with it, both the default
            coordinate system type and code must be specified for a coordinate system to be assigned to
            the returning JSON. The original geometry is left unaltered.
        """
        ref_sys = geom.GetSpatialReference()
        geom_json = json.loads(geom.ExportToJson())
        if not ref_sys:
            if alt_coord_type and alt_coord_code:
                geom_json['crs'] = {'type': str(alt_coord_type), 'properties': {'code': str(alt_coord_code)}}
        else:
            geom_json['crs'] = {
                'type': ref_sys.GetAttrValue("AUTHORITY", 0),
                'properties': {
                    'code': ref_sys.GetAttrValue("AUTHORITY", 1)
                }
            }

        return json.dumps(geom_json)

    @staticmethod
    def convert_geometry(geometry: ogr.Geometry, new_spatialreference: osr.SpatialReference) -> ogr.Geometry:
        """Converts the geometry to the new spatial reference if possible

        geometry - The geometry to transform
        new_spatialreference - The spatial reference to change to

        Returns:
            The transformed geometry or the original geometry. If either the
            new Spatial Reference parameter is None, or the geometry doesn't
            have a spatial reference, then the original geometry is returned.
        """
        if not new_spatialreference or not geometry:
            return geometry

        return_geometry = geometry
        try:
            geom_sr = geometry.GetSpatialReference()
            if geom_sr and not new_spatialreference.IsSame(geom_sr):
                transform = osr.CreateCoordinateTransformation(geom_sr, new_spatialreference)
                new_geom = geometry.Clone()
                if new_geom:
                    new_geom.Transform(transform)
                    return_geometry = new_geom
        except Exception as ex:
            logging.warning("Exception caught while transforming geometries: %s", str(ex))
            logging.warning("    Returning original geometry")

        return return_geometry

    @staticmethod
    def find_plots_intersect_boundingbox(bounding_box: str, all_plots: dict, fullmac: bool = True) -> dict:
        """Take a list of plots from BETY and return only those overlapping bounding box.
        Arguments:
            bounding_box: the JSON of the bounding box
            all_plots: the dictionary of all available plots
            fullmac: only include full plots (omit KSU, omit E W partial plots)
        Return:
            A dictionary of all intersecting plots
        """
        bbox_poly = ogr.CreateGeometryFromJson(str(bounding_box))
        bb_sr = bbox_poly.GetSpatialReference()
        intersecting_plots = {}
        logging.debug("HACK: Bounding box %s %s", str(bb_sr), str(bbox_poly))

        for plotname in all_plots:
            if fullmac and (plotname.find("KSU") > -1 or plotname.endswith(" E") or plotname.endswith(" W")):
                continue

            bounds = all_plots[plotname]

            yaml_bounds = yaml.safe_load(bounds)
            current_poly = ogr.CreateGeometryFromJson(json.dumps(yaml_bounds))

            # Check for a need to convert coordinate systems
            check_poly = current_poly
            if bb_sr:
                poly_sr = current_poly.GetSpatialReference()
                if poly_sr and not bb_sr.IsSame(poly_sr):
                    # We need to convert to the same coordinate system before an intersection
                    check_poly = __internal__.convert_geometry(current_poly, bb_sr)

            logging.debug("HACK: Intersection with %s", str(check_poly))
            intersection_with_bounding_box = bbox_poly.Intersection(check_poly)

            if intersection_with_bounding_box is not None:
                intersection = json.loads(intersection_with_bounding_box.ExportToJson())
                if 'coordinates' in intersection and len(intersection['coordinates']) > 0:
                    intersecting_plots[plotname] = bounds

        return intersecting_plots

    @staticmethod
    def image_get_geobounds(filename: str) -> list:
        """Uses gdal functionality to retrieve recilinear boundaries from the file

        Args:
            filename(str): path of the file to get the boundaries from

        Returns:
            The upper-left and calculated lower-right boundaries of the image in a list upon success.
            The values are returned in following order: min_y, max_y, min_x, max_x. A list of numpy.nan
            in each position is returned if the boundaries can't be determined
        """
        try:
            src = gdal.Open(filename)
            ulx, xres, _, uly, _, yres = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)

            min_y = min(uly, lry)
            max_y = max(uly, lry)
            min_x = min(ulx, lrx)
            max_x = max(ulx, lrx)

            return [min_y, max_y, min_x, max_x]
        except Exception as ex:
            logging.info("[image_get_geobounds] Exception caught: %s", str(ex))
            if logging.getLogger().level == logging.DEBUG:
                logging.exception("[image_get_geobounds] Exception")

        return [np.nan, np.nan, np.nan, np.nan]

    @staticmethod
    def clip_raster(rast_path: str, bounds: tuple, out_path: str = None, nodata: int = -9999, compress: bool = False) ->\
            Optional[np.ndarray]:
        """Clip raster to polygon.
        Arguments:
            rast_path: path to raster file
            bounds: (min_y, max_y, min_x, max_x)
            out_path: if provided, where to save as output file
            nodata: the no data value
            compress: set to True to compress the image and False to leave it uncompressed
        Returns:
            The array of clipped pixels
        Notes:
            Oddly, the "features path" can be either a filename
            OR a geojson string. GDAL seems to figure it out and do
            the right thing.
            From http://karthur.org/2015/clipping-rasters-in-python.html
        """
        if not out_path:
            out_path = "temp.tif"

        # Clip raster to GDAL and read it to numpy array
        coords = "%s %s %s %s" % (bounds[2], bounds[1], bounds[3], bounds[0])
        if compress:
            cmd = 'gdal_translate -projwin %s "%s" "%s"' % (coords, rast_path, out_path)
        else:
            cmd = 'gdal_translate -co COMPRESS=LZW -projwin %s "%s" "%s"' % (coords, rast_path, out_path)
        subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        out_px = np.array(gdal.Open(out_path).ReadAsArray())

        if np.count_nonzero(out_px) > 0:
            if out_path == "temp.tif":
                os.remove(out_path)
            return out_px

        os.remove(out_path)
        return None

    @staticmethod
    def get_epsg(filename: str) -> Optional[str]:
        """Returns the EPSG of the georeferenced image file
        Args:
            filename(str): path of the file to retrieve the EPSG code from
        Return:
            Returns the found EPSG code, or None if it's not found or an error ocurred
        """
        logger = logging.getLogger(__name__)

        try:
            src = gdal.Open(filename)

            proj = osr.SpatialReference(wkt=src.GetProjection())

            return proj.GetAttrValue('AUTHORITY', 1)
        except Exception as ex:
            logger.warning("[get_epsg] Exception caught: %s", str(ex))
            if logging.getLogger().level == logging.DEBUG:
                logging.exception("[get_epsg] Exception")

        return None

    @staticmethod
    def get_las_epsg_from_header(header: liblas.header.Header) -> str:
        """Returns the found EPSG code from the LAS header
        Arguments:
            header: the loaded LAS header to find the SRID in
        Return:
            Returns the SRID as a string if found, None is returned otherwise
        """
        epsg = None
        search_terms_ordered = ['DATUM', 'AUTHORITY', '"EPSG"', ',']
        try:
            # Get the WKT from the header, find the DATUM, then finally the EPSG code
            srs = header.get_srs()
            wkt = srs.get_wkt().decode('UTF-8')
            idx = -1
            for term in search_terms_ordered:
                idx = wkt.find(term)
                if idx < 0:
                    break
            if idx >= 0:
                epsg = re.search(r'\d+', wkt[idx:])[0]
        except Exception as ex:
            logging.debug("Unable to find EPSG in LAS file header")
            logging.debug("    exception caught: %s", str(ex))

        return epsg

    @staticmethod
    def get_las_extents(file_path: str, default_epsg: int = None) -> Optional[str]:
        """Calculate the extent of the given las file and return as GeoJSON.
        Arguments:
            file_path: path to the file from which to load the bounds
            default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
        Return:
            Returns the JSON representing the image boundary, or None if the
            bounds could not be loaded
        Notes:
            If a file doesn't have a coordinate system and a default epsg is specified, the
            return JSON will use the default_epsg.
            If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
            of the image is not returned (None) and a warning is logged.
        """
        # Get the bounds and the EPSG code
        las_info = liblas.file.File(file_path, mode='r')
        min_bound = las_info.header.min
        max_bound = las_info.header.max
        epsg = __internal__.get_las_epsg_from_header(las_info.header)
        if epsg is None:
            if default_epsg is not None:
                epsg = default_epsg
            else:
                logging.warning("Unable to find EPSG and not default is specified for file '%s'", file_path)
                return None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(min_bound[1], min_bound[0])  # Upper left
        ring.AddPoint(min_bound[1], max_bound[0])  # Upper right
        ring.AddPoint(max_bound[1], max_bound[0])  # lower right
        ring.AddPoint(max_bound[1], min_bound[0])  # lower left
        ring.AddPoint(min_bound[1], min_bound[0])  # Closing the polygon

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        ref_sys = osr.SpatialReference()
        if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
            poly.AssignSpatialReference(ref_sys)
            return __internal__.geometry_to_geojson(poly)

        logging.error("Failed to import EPSG %s for las file %s", str(epsg), file_path)
        return None

    @staticmethod
    def clip_las(las_path: str, clip_tuple: tuple, out_path: str) -> None:
        """Clip LAS file to polygon.
        Arguments:
          las_path: path to point cloud file
          clip_tuple: tuple containing (minX, maxX, minY, maxY) of clip bounds
          out_path: output file to write
        Notes:
            The clip_tuple is assumed to be in the correct coordinate system for the point cloud file
        """
        bounds_str = "([%s, %s], [%s, %s])" % (clip_tuple[0], clip_tuple[1], clip_tuple[2], clip_tuple[3])

        pdal_dtm = out_path.replace(".las", "_dtm.json")
        with open(pdal_dtm, 'w') as dtm:
            dtm_data = """{
                "pipeline": [
                    "%s",
                    {
                        "type": "filters.crop",
                        "bounds": "%s"
                    },
                    {
                        "type": "writers.las",
                        "filename": "%s"
                    }
                ]
            }""" % (las_path, bounds_str, out_path)
            logging.debug("Writing dtm file contents: %s", str(dtm_data))
            dtm.write(dtm_data)

        cmd = 'pdal pipeline "%s"' % pdal_dtm
        logging.debug("Running pipeline command: %s", cmd)
        subprocess.call([cmd], shell=True)
        os.remove(pdal_dtm)

    @staticmethod
    def get_image_bounds_json(file_path: str, default_epsg: int = None) -> Optional[str]:
        """Loads the boundaries of the image file and returns the GeoJSON
           representing the bounds (including EPSG code)
        Arguments:
            file_path: path to the file from which to load the bounds
            default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
        Return:
            Returns the JSON representing the image boundary, or None if the
            bounds could not be loaded
        Notes:
            If a file doesn't have a coordinate system and a default epsg is specified, the
            return JSON will use the default_epsg.
            If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
            of the image is not returned (None) and a warning is logged.
        """
        # Get the bounds (if they exist)
        bounds = __internal__.image_get_geobounds(file_path)
        if bounds[0] == np.nan:
            return None

        epsg = __internal__.get_epsg(file_path)
        if epsg is None:
            if default_epsg:
                epsg = default_epsg
            else:
                logging.warning("Files does not have a coordinate system defined and no default was specified: '%s'",
                                file_path)
                return None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(bounds[2], bounds[1])  # Upper left
        ring.AddPoint(bounds[3], bounds[1])  # Upper right
        ring.AddPoint(bounds[3], bounds[0])  # lower right
        ring.AddPoint(bounds[2], bounds[0])  # lower left
        ring.AddPoint(bounds[2], bounds[1])  # Closing the polygon

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        ref_sys = osr.SpatialReference()
        if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
            poly.AssignSpatialReference(ref_sys)
            return __internal__.geometry_to_geojson(poly)

        logging.error("Failed to import EPSG %s for image file %s", str(epsg), file_path)
        return None

    @staticmethod
    def get_files_to_process(file_list: list, sensor: str, default_epsg: int = None) -> dict:
        """Returns a dictionary of georeferenced files to process
        Arguments:
            file_list: the list of file paths to process
            sensor: the name of the sensor associated with the files
            default_epsg: the default EPSG value to use if a file is missing one
        Return:
            Returns a dictionary with the file names as keys. Each key's value is another dictionary containing
            the file path, file bounds (as GeoJSON), and the sensor name
        """
        files_to_process = {}
        for one_file in file_list:
            filename = os.path.basename(one_file)
            if filename in files_to_process:
                continue
            if not os.path.exists(one_file):
                logging.warning("Skipping file that does not exist: '%s'", one_file)
                continue

            if one_file.endswith('.tif'):
                files_to_process[filename] = {
                    'path': one_file,
                    'bounds': __internal__.get_image_bounds_json(one_file, default_epsg),
                    'sensor_name': sensor
                }
            elif one_file.endswith(".las"):
                files_to_process[filename] = {
                    'path': one_file,
                    'bounds': __internal__.get_las_extents(one_file, default_epsg),
                    'sensor_name': sensor
                }
        return files_to_process

    @staticmethod
    def get_spatial_reference_from_json(geojson: str) -> Optional[osr.SpatialReference]:
        """Returns the spatial reference embedded in the geojson.
        Args:
            geojson(str): the geojson to get the spatial reference from
        Return:
            The osr.SpatialReference that represents the geographic coordinate system
            in the geojson. None is returned if a spatial reference isn't found
        """
        yaml_geom = yaml.safe_load(geojson)
        current_geom = ogr.CreateGeometryFromJson(json.dumps(yaml_geom))

        if current_geom:
            return current_geom.GetSpatialReference()
        return None

    @staticmethod
    def geojson_to_tuples(bounding_box: str) -> tuple:
        """Returns the bounds of the shape
        Arguments:
            bounding_box: the JSON of the geometry
        Return:
            A tuple containing the bounds in (min Y, max Y, min X, max X) order
        """
        yaml_geom = yaml.safe_load(bounding_box)
        current_geom = ogr.CreateGeometryFromJson(json.dumps(yaml_geom))
        current_env = current_geom.GetEnvelope()

        return current_env[2], current_env[3], current_env[0], current_env[1]

    @staticmethod
    def calculate_overlap_percent(check_bounds: str, bounding_box: str) -> float:
        """Calculates and returns the percentage overlap between the two boundaries.
           The calculation determines the overlap shape between the two parameters and
           then calculates the percentage by dividing the overlap area by the bounding
           box area, and returns that value.
        Args:
            check_bounds: GeoJSON of boundary to check
            bounding_box: GeoJSON of boundary to check against
        Return:
            The calculated overlap percent (0.0 - 1.0) or 0.0 if there is no overlap.
            If an exception is detected, a warning message is logged and 0.0 is returned.
        """
        try:
            check_poly = ogr.CreateGeometryFromJson(str(check_bounds))
            bbox_poly = ogr.CreateGeometryFromJson(str(bounding_box))

            if check_poly and bbox_poly:
                intersection = bbox_poly.Intersection(check_poly)
                if intersection:
                    return intersection.Area() / check_poly.Area()
        except Exception as ex:
            logging.warning("Exception caught while calculating shape overlap: %s", str(ex))

        return 0.0

    @staticmethod
    def clip_raster_intersection(file_path: str, file_bounds: str, plot_bounds: str, out_file: str) -> Optional[int]:
        """Clips the raster to the intersection of the file bounds and plot bounds
        Arguments:
            file_path: the path to the source file
            file_bounds: the geometric boundary of the source file as JSON
            plot_bounds: the geometric boundary of the plot to clip to as JSON
            out_file: the path to store the clipped image
        Return:
            The number of pixels in the new image, or None if no pixels were saved
        Notes:
            Assumes the boundaries are in the same coordinate system
        Exceptions:
            Raises RuntimeError if the polygons are invalid
        """
        logging.debug("Clip to intersect of plot boundary: File: '%s' '%s' Plot: '%s'", file_path, str(file_bounds), str(plot_bounds))
        try:
            file_poly = ogr.CreateGeometryFromJson(str(file_bounds))
            plot_poly = ogr.CreateGeometryFromJson(str(plot_bounds))

            if not file_poly or not plot_poly:
                logging.error("Invalid polygon specified for clip_raster_intersection: File: '%s' plot: '%s'",
                              str(file_bounds), str(plot_bounds))
                raise RuntimeError("One or more invalid polygons specified when clipping raster")

            intersection = file_poly.Intersection(plot_poly)
            if not intersection or not intersection.Area():
                logging.info("File does not intersect plot boundary: %s", file_path)
                return None

            # Make sure we pass a multipolygon down to the tuple converter
            if intersection.GetGeometryName().startswith('MULTI'):
                multi_polygon = intersection
            else:
                multi_polygon = ogr.Geometry(ogr.wkbMultiPolygon)
                multi_polygon.AddGeometry(intersection)

            # Proceed to clip to the intersection
            tuples = __internal__.geojson_to_tuples(__internal__.geometry_to_geojson(multi_polygon))
            return __internal__.clip_raster(file_path, tuples, out_path=out_file, compress=True)

        except Exception as ex:
            logging.exception("Exception caught while clipping image to plot intersection")
            raise ex

    @staticmethod
    def cleanup_request_md(source_md: dict) -> dict:
        """Makes a copy of the source metadata and cleans it up for use as plot-level information
        Arguments:
            source_md: the source metadata to clone and clean up
        Returns:
            returns the cleaned up metadata
        """
        if not source_md:
            return {}

        new_md = copy.deepcopy(source_md)
        new_md.pop('list_files', None)
        new_md.pop('context_md', None)
        new_md.pop('working_folder', None)

        return new_md

    @staticmethod
    def check_already_merged(merged_file: str, source_file: str) -> bool:
        """Checks if the source file is listed in the merged file
        Arguments:
            merged_file: path to the merged file to check
            source_file: the name of the source file to look for
        Return:
            Returns True if the source file name is found in the contents of the merged file and False otherwise
        """
        already_merged = False
        if os.path.exists(merged_file):
            # Check if contents
            with open(merged_file, 'r') as contents:
                for entry in contents.readlines():
                    if entry.strip() == source_file:
                        already_merged = True
                        break
        return already_merged

    @staticmethod
    def prepare_container_md(plot_name: str, plot_md: dict, sensor: str, source_file: str, result_files: list) -> dict:
        """Prepares the metadata for a single container
        Arguments:
            plot_name: the name of the container
            plot_md: the metadata associated with this container
            sensor: the name of the related sensor
            source_file: the name of the source file
            result_files: list of files to add to container metadata
        Return:
            The formatted metadata
        Notes:
            The files in result_files are checked for existence before being added to the metadata
        """
        cur_md = {
            'name': plot_name,
            'metadata': {
                'replace': True,
                'data': plot_md
            },
            'file': []
        }
        for one_file in result_files:
            if os.path.exists(one_file):
                cur_md['file'].append({
                    'path': one_file,
                    'key': sensor,
                    'metadata': {
                        'source': source_file,
                        'transformer': configuration.TRANSFORMER_NAME,
                        'version': configuration.TRANSFORMER_VERSION,
                        'timestamp': datetime.datetime.utcnow().isoformat(),
                        'plot_name': plot_name
                    }
                })
        return cur_md

    @staticmethod
    def merge_container_md(dest_md: list, new_md: dict) -> list:
        """Merges container level metadata ensuring there aren't any plot-level
           duplicates or duplicate file entries for a plot entry
        Arguments:
            dest_md: the list of current metadata to merge into
            new_md: the new metadata to merge
        Return:
            Returns a new list of metadata with the new metadata merged into it
        """
        # Return something meaningful if we have missing or empty dict
        if not dest_md:
            if new_md:
                return [new_md]
            return []

        # Copy the metadata and look for a match
        match_idx = -1
        md_len = len(dest_md)
        for idx in range(0, md_len):
            if dest_md[idx]['name'] == new_md['name']:
                match_idx = idx
                break

        # If no match found, add and return
        if match_idx == -1:
            dest_md.append(new_md)
            return dest_md

        # Merge the metadata
        working_md = dest_md[match_idx]
        if 'files' in new_md:
            if 'files' in working_md:
                # Only add files that aren't included in the destination metadata already
                for one_file in new_md['files']:
                    file_match_found = False
                    for match_file in working_md['files']:
                        if one_file['path'] == match_file['path']:
                            file_match_found = True
                            break
                    if not file_match_found:
                        dest_md[match_idx]['files'].append(one_file)
            else:
                # Target metadata doesn't have a 'files' entry
                dest_md[match_idx]['files'] = new_md['files']

        return dest_md


def add_parameters(parser: argparse.ArgumentParser) -> None:
    """Adds parameters
    Arguments:
        parser: instance of argparse
    """
    parser.add_argument('--epsg', type=int, nargs='?',
                        help='default epsg code to use if a file doesn\'t have a coordinate system')
    parser.add_argument('--full_plot_fill', action='store_true',
                        help='clipped images will be color filled to match the plot dimensions (outside the original image boundaries)')
    parser.add_argument('--shapefile', type=str, help='the path of the shapefile to use for plot boundaries')
    parser.add_argument('--shapefile_plot_column', type=str,
                        help='the name of the column in the shapefile containing plot names (defaults to "plot_" + row number)')
    parser.add_argument('sensor', type=str, help='the name of the sensor associated with the source files')


def perform_process(transformer: transformer_class.Transformer, check_md: dict, transformer_md: dict, full_md: list) -> dict:
    """Performs the processing of the data
    Arguments:
        transformer: instance of transformer class
        check_md: metadata associated with this request
        transformer_md: metadata associated with this transformer
        full_md: the full set of metadata
    Return:
        Returns a dictionary with the results of processing
    """
    # pylint: disable=unused-argument
    # loop through the available files and clip data into plot-level files
    processed_files = 0
    processed_plots = 0
    start_timestamp = datetime.datetime.now()
    file_list = check_md['list_files']()
    files_to_process = __internal__.get_files_to_process(file_list, transformer.args.sensor, transformer.args.epsg)
    logging.info("Found %s files to process", str(len(files_to_process)))

    container_md = []
    if files_to_process:
        # Get all the possible plots
        datestamp = check_md['timestamp'][0:10]
        if transformer.args.shapefile:
            all_plots = __internal__.load_shapefile(transformer.args.shapefile, transformer.args.shapefile_plot_column)
        else:
            all_plots = __internal__.get_site_boundaries(datestamp, city='Maricopa')
            logging.debug("Have %s plots for site", len(all_plots))

        for filename in files_to_process:
            processed_files += 1
            file_path = files_to_process[filename]['path']
            file_bounds = files_to_process[filename]['bounds']
            sensor = files_to_process[filename]['sensor_name']
            logging.debug("File bounds: %s", str(file_bounds))

            overlap_plots = __internal__.find_plots_intersect_boundingbox(file_bounds, all_plots, fullmac=True)
            logging.info("Have %s plots intersecting file '%s'", str(len(overlap_plots)), filename)

            file_spatial_ref = __internal__.get_spatial_reference_from_json(file_bounds)
            for plot_name in overlap_plots:
                processed_plots += 1
                plot_bounds = __internal__.convert_json_geometry(overlap_plots[plot_name], file_spatial_ref)
                logging.debug("Clipping out plot '%s': %s", str(plot_name), str(plot_bounds))
                if __internal__.calculate_overlap_percent(plot_bounds, file_bounds) < 0.10:
                    logging.info("Skipping plot with too small overlap: %s", plot_name)
                    continue
                tuples = __internal__.geojson_to_tuples(plot_bounds)

                plot_md = __internal__.cleanup_request_md(check_md)
                plot_md['plot_name'] = plot_name

                if filename.endswith('.tif'):
                    # If file is a geoTIFF, simply clip it
                    out_path = os.path.join(check_md['working_folder'], plot_name)
                    out_file = os.path.join(out_path, filename)
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)

                    if not transformer.args.full_plot_fill:
                        __internal__.clip_raster_intersection(file_path, file_bounds, plot_bounds, out_file)
                    else:
                        logging.info("Clipping image to plot boundary with fill")
                        __internal__.clip_raster(file_path, tuples, out_path=out_file, compress=True)

                    cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file])
                    container_md = __internal__.merge_container_md(container_md, cur_md)

                elif filename.endswith('.las'):
                    out_path = os.path.join(check_md['working_folder'], plot_name)
                    out_file = os.path.join(out_path, filename)
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)

                    __internal__.clip_las(file_path, tuples, out_path=out_file)

                    cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file])
                    container_md = __internal__.merge_container_md(container_md, cur_md)

    return {
        'code': 0,
        'container': container_md,
        configuration.TRANSFORMER_NAME:
        {
            'utc_timestamp': datetime.datetime.utcnow().isoformat(),
            'processing_time': str(datetime.datetime.now() - start_timestamp),
            'total_file_count': len(file_list),
            'processed_file_count': processed_files,
            'total_plots_processed': processed_plots,
            'sensor': transformer.args.sensor
        }
    }
