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
import netCDF4

import yaml
import liblas

from osgeo import ogr
import numpy as np
import osr

from terrautils.imagefile import image_get_geobounds, get_epsg
from terrautils.spatial import find_plots_intersect_boundingbox, geojson_to_tuples_betydb, clip_raster, \
     geometry_to_geojson, convert_json_geometry
from terrautils.betydb import get_site_boundaries
import terrautils.lemnatec

import configuration
import transformer_class

terrautils.lemnatec.SENSOR_METADATA_CACHE = os.path.dirname(os.path.realpath(__file__))

# The default EPSG code expected for netCDF files
DEFAULT_NETCDF_EPSG_CODE = 4326


class __internal__():
    """Class for internal use only functions
    """

    def __init__(self):
        """Initializes class instance
        """

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
            return geometry_to_geojson(poly)

        logging.error("Failed to import EPSG %s for las file %s", str(epsg), file_path)
        return None

    @staticmethod
    def get_netcdf_extents(file_path: str) -> Optional[str]:
        """Calculate the extent of the given netCDF file and return as GeoJSON.
        Arguments:
            file_path: path to the file from which to load the bounds
        Return:
            Returns the JSON representing the netCDF boundary, or None if the
            bounds could not be loaded
        """
        epsg = DEFAULT_NETCDF_EPSG_CODE
        try:
            with netCDF4.Dataset(file_path) as src:
                # Get the variables we need (assumes a 90 degree rotation from North-up)
                y_size = src.dimensions['y'].size
                x_size = src.dimensions['x'].size
                lat_min = src.variables['latitude'][0:1].data[0]
                lon_min = src.variables['longitude'][0:1].data[0]
                lat_max = src.variables['latitude'][x_size - 1:x_size].data[0]
                lon_max = src.variables['longitude'][y_size - 1:y_size].data[0]

                if lat_min > lat_max:
                    lat_min, lat_max = lat_max, lat_min
                if lon_min > lon_max:
                    lon_min, lon_max = lon_max, lon_min

                ring = ogr.Geometry(ogr.wkbLinearRing)
                ring.AddPoint(lon_min, lat_max)  # Upper left
                ring.AddPoint(lon_max, lat_max)  # Upper right
                ring.AddPoint(lon_max, lat_min)  # lower right
                ring.AddPoint(lon_min, lat_min)  # lower left
                ring.AddPoint(lon_min, lat_max)  # Closing the polygon

                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)

                ref_sys = osr.SpatialReference()
                if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
                    poly.AssignSpatialReference(ref_sys)
                    return geometry_to_geojson(poly)

                logging.error("Failed to import EPSG %s for las file %s", str(epsg), file_path)
        except Exception as ex:
            logging.exception("Exception caught while determining netCDF file extents %s", file_path)
            raise ex

        return None

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
        bounds = image_get_geobounds(file_path)
        if bounds[0] == np.nan:
            return None

        epsg = get_epsg(file_path)
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
            return geometry_to_geojson(poly)

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
            elif one_file.endswith(".nc"):
                files_to_process[filename] = {
                    'path': one_file,
                    'bounds': __internal__.get_netcdf_extents(one_file),
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
    def netcdf_tuple_to_indexes(netcdf_path: str, clip_tuple: tuple) -> list:
        """Converts the tuple lat-lon values to Y and X indexes
        Arguments:
          netcdf_path: path to netCDF file
          clip_tuple: tuple containing (minX, maxX, minY, maxY) to convert
        Return:
            Returns a tuple containing (minYindex, maxYimdex, minXindex, maxXindex)
        Exceptions:
            Raises RuntimeError if the tuple is not able to be mapped to X, Y indexes
        Notes:
            The X, Y axises of the file are aligned east-west, north-south (x is oriented along latitude)
        """
        closest_offsets = []
        with netCDF4.Dataset(netcdf_path) as src:
            # Get the variables we need (assumes a 90 degree rotation from North-up)
            y_size = float(src.dimensions['y'].size)
            x_size = float(src.dimensions['x'].size)
            latitudes = src.variables['latitude']
            longitudes = src.variables['longitude']
            lat_min = latitudes[0:1].data[0]
            lon_min = longitudes[0:1].data[0]
            lat_max = latitudes[x_size - 1:x_size].data[0]
            lon_max = longitudes[y_size - 1:y_size].data[0]

            # Get approximate X and Y of tuple values - assumes constant spacing of sensor points
            # Y is assumed to be aligned along longitude
            lat_diff = abs(lat_max - lat_min)
            lon_diff = abs(lon_max - lon_min)
            min_y = int(y_size * (abs(clip_tuple[0] - lon_min) / lon_diff))
            max_y = int(y_size * (abs(clip_tuple[1] - lon_max) / lon_diff))
            min_x = int(x_size * (abs(clip_tuple[2] - lat_min) / lat_diff))
            max_x = int(x_size * (abs(clip_tuple[3] - lat_max) / lat_diff))
            logging.debug("Approximate indexes: %s %s %s %s", str(min_y), str(max_y), str(min_x), str(max_x))

            # Find the closest offsets to the plot boundary
            def find_closest_index(value: float, value_index: int, variables: netCDF4.Variable) -> int:
                """Finds the closes indexes for the values in the list
                Arguments:
                    value: the value to look up
                    value_index: the index of value currently deemed closest
                    variables: the variables to look into - assumed to be in order
                Return:
                    Returns an index corresponding to the nearest match to the values. If a value can't
                    be matched a None is returned
                """
                variables_size = len(variables)
                closest = None
                min_diff = 9999999999.99
                for index_offset in range(-2, 3):
                    test_index = value_index + index_offset
                    if test_index < 0 or test_index >= variables_size:
                        continue
                    test_variable = variables[test_index:test_index + 1].data[0]
                    cur_diff = abs(test_variable - value)
                    if cur_diff < min_diff:
                        closest = test_index
                        min_diff = cur_diff

                return closest

            min_y_nearest = find_closest_index(lon_min, min_y, longitudes)
            if None in min_y_nearest:
                raise RuntimeError("Unable to find index for minimum Y location value")
            max_y_nearest = find_closest_index(lon_max, max_y, longitudes)
            if None in max_y_nearest:
                raise RuntimeError("Unable to find index for maximum Y location value")
            min_x_nearest = find_closest_index(lat_min, min_x, latitudes)
            if None in min_x_nearest:
                raise RuntimeError("Unable to find index for minimum X location value")
            max_x_nearest = find_closest_index(lat_max, max_x, latitudes)
            if None in max_x_nearest:
                raise RuntimeError("Unable to find index for maximum X location value")

            if min_y_nearest > max_y_nearest:
                min_y_nearest, max_y_nearest = max_y_nearest, min_y_nearest
            if min_x_nearest > max_x_nearest:
                min_x_nearest, max_x_nearest = max_x_nearest, min_x_nearest

            closest_offsets.append(min_y_nearest)
            closest_offsets.append(max_y_nearest)
            closest_offsets.append(min_x_nearest)
            closest_offsets.append(max_x_nearest)

        logging.debug("Closest indexes: %s", str(closest_offsets))
        return closest_offsets

    @staticmethod
    def clip_netcdf(netcdf_path: str, clip_tuple: tuple, out_path: str) -> None:
        """Clip netCDF file to polygon.
        Arguments:
          netcdf_path: path to netCDF file
          clip_tuple: tuple containing (minX, maxX, minY, maxY) of clip bounds
          out_path: output file to write
        Notes:
            The clip_tuple is assumed to be in the correct coordinate system of the netCDF file
        """
        logging.debug("Clipping netCDF file to plot tuple: %s", str(clip_tuple))

        # Map the clip tuple to X and Y positions
        closest_offsets = __internal__.netcdf_tuple_to_indexes(netcdf_path, clip_tuple)

        # Make the NCO calls to clip
        cmd = 'ncks'
        command_params = ['-O',
                          '-d',
                          'x,{},{}'.format(closest_offsets[2], closest_offsets[3]),
                          '-d y,{},{}'.format(closest_offsets[0], closest_offsets[1]),
                          netcdf_path,
                          out_path]
        logging.debug("Running netCDF clipping command: '%s' parameters: %s", cmd, str(command_params))
        run_result = subprocess.run([cmd, *command_params], check=True)
        logging.debug("Command result: %s", run_result.returncode)
        logging.debug("        output: %s", str(run_result.stdout))
        logging.debug("         error: %s", str(run_result.stderr))

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

    def convert_geometry(geometry, new_spatialreference):
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
            logging.warning("Exception caught while transforming geometries: " + str(ex))
            logging.warning("    Returning original geometry")

        return return_geometry

    def find_plots_intersect_boundingbox(bounding_box, all_plots, fullmac=True):
        """Take a list of plots from BETY and return only those overlapping bounding box.

        fullmac -- only include full plots (omit KSU, omit E W partial plots)

        """
        bbox_poly = ogr.CreateGeometryFromJson(str(bounding_box))
        bb_sr = bbox_poly.GetSpatialReference()
        intersecting_plots = dict()

        for plotname in all_plots:
            if fullmac and (plotname.find("KSU") > -1 or plotname.endswith(" E") or plotname.endswith(" W")):
                continue

            bounds = all_plots[plotname]
            if 'Season 4 Range 37 Column 15' in plotname:
                logging.debug("Checking plot %s: %s", plotname, str(bounds))

            yaml_bounds = yaml.safe_load(bounds)
            current_poly = ogr.CreateGeometryFromJson(json.dumps(yaml_bounds))

            # Check for a need to convert coordinate systems
            check_poly = current_poly
            if bb_sr:
                poly_sr = current_poly.GetSpatialReference()
                if poly_sr and not bb_sr.IsSame(poly_sr):
                    # We need to convert to the same coordinate system before an intersection
                    logging.debug("    Converting to new coordinate system")
                    check_poly = __internal__.convert_geometry(current_poly, bb_sr)

            intersection_with_bounding_box = bbox_poly.Intersection(check_poly)

            if intersection_with_bounding_box is not None:
                intersection = json.loads(intersection_with_bounding_box.ExportToJson())
                if 'coordinates' in intersection and len(intersection['coordinates']) > 0:
                    logging.debug("  HAVE INTERSECTION")
                    intersecting_plots[plotname] = bounds

        return intersecting_plots


def add_parameters(parser: argparse.ArgumentParser) -> None:
    """Adds parameters
    Arguments:
        parser: instance of argparse
    """
    parser.add_argument('--epsg', type=int, nargs='?',
                        help='default epsg code to use if a files doesn\'t have a coordinate system')
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

    # Get all the possible plots
    datestamp = check_md['timestamp'][0:10]
    logging.debug("Using datestamp %s (%s)", datestamp, check_md['timestamp'])
    all_plots = get_site_boundaries(datestamp, city='Maricopa')
    logging.debug("Have %s plots for site", len(all_plots))

    container_md = []
    for filename in files_to_process:
        processed_files += 1
        file_path = files_to_process[filename]['path']
        file_bounds = files_to_process[filename]['bounds']
        sensor = files_to_process[filename]['sensor_name']
        logging.debug("File bounds: %s", str(file_bounds))
        logging.debug("First plot: %s", str(all_plots)[0:1000])

        overlap_plots = __internal__.find_plots_intersect_boundingbox(file_bounds, all_plots, fullmac=True)
        logging.info("Have %s plots intersecting file '%s'", str(len(overlap_plots)), filename)

        file_spatial_ref = __internal__.get_spatial_reference_from_json(file_bounds)
        for plot_name in overlap_plots:
            processed_plots += 1
            plot_bounds = convert_json_geometry(overlap_plots[plot_name], file_spatial_ref)
            logging.debug("Clipping out plot '%s': %s", str(plot_name), str(plot_bounds))
            if __internal__.calculate_overlap_percent(plot_bounds, file_bounds) < 0.10:
                logging.info("Skipping plot with too small overlap: %s", plot_name)
                continue
            tuples = geojson_to_tuples_betydb(yaml.safe_load(plot_bounds))

            plot_md = __internal__.cleanup_request_md(check_md)
            plot_md['plot_name'] = plot_name

            # Prepare for clipping. We use a lambda function to determine if we're actually clipping a file
            # and to perform the clipping allowing us to simplify the code
            clip_func = None
            out_path = os.path.join(check_md['working_folder'], plot_name)
            out_file = os.path.join(out_path, filename)
            # Silence the following check since we run the lambda as soon as possible after it's defined
            # pylint: disable=cell-var-from-loop
            if filename.endswith('.tif'):
                # If file is a geoTIFF, simply clip it
                clip_func = lambda: clip_raster(file_path, tuples, out_path=out_file, compress=True)
            elif filename.endswith('.las'):
                clip_func = lambda: __internal__.clip_las(file_path, tuples, out_file)
            elif filename.endswith('.nc'):
                clip_func = lambda: __internal__.clip_netcdf(file_path, tuples, out_file)

            if clip_func:
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                clip_func()

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
