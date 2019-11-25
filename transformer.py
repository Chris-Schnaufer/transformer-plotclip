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
    def get_las_extents(file_path: str, default_epsg: int = None) -> str:
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
    def clip_las(las_path: str, clip_tuple: tuple, out_path: str, merged_path: str = None) -> None:
        """Clip LAS file to polygon.
        Arguments:
          las_path: path to pointcloud file
          clip_tuple: tuple containing (minX, maxX, minY, maxY) of clip bounds
          out_path: output file to write
          merge_path: optional path to write merged data to
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

        if merged_path:
            if os.path.isfile(merged_path):
                cmd = 'pdal merge "%s" "%s" "%s"' % (out_path, merged_path, merged_path)
                logging.debug("Running merge command: %s", cmd)
                subprocess.call([cmd], shell=True)
            else:
                os.rename(out_path, merged_path)

    @staticmethod
    def get_image_bounds_json(file_path: str, default_epsg: int = None) -> str:
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
        return files_to_process

    @staticmethod
    def get_spatial_reference_from_json(geojson):
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
    def calculate_overlap_percent(check_bounds, bounding_box):
        """Calculates and returns the percentage overlap between the two boundaries.
           The calculation determines the overlap shape between the two parameters and
           then calculates the percentage by dividing the overlap area by the bounding
           box area, and returns that value.
        Args:
            check_bounds(str): GeoJSON of boundary to check
            bounding_box(str): GeoJSON of boundary to check against
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
            return {}

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
                        help='default epsg code to use if a files doesn\'t have a coordinate system')
    parser.add_argument('sensor', type=str, help='the name of the sensor associated with the source files')


def perform_process(transformer: transformer_class.Transformer, check_md: dict, transformer_md: dict, full_md: dict) -> dict:
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
    all_plots = get_site_boundaries(datestamp, city='Maricopa')
    logging.debug("Have %s plots for site", len(all_plots))

    container_md = []
    for filename in files_to_process:
        processed_files += 1
        file_path = files_to_process[filename]['path']
        file_bounds = files_to_process[filename]['bounds']
        sensor = files_to_process[filename]['sensor_name']

        overlap_plots = find_plots_intersect_boundingbox(file_bounds, all_plots, fullmac=True)
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

            if filename.endswith('.tif'):
                # If file is a geoTIFF, simply clip it
                out_path = os.path.join(check_md['working_folder'], plot_name)
                out_file = os.path.join(out_path, filename)
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                clip_raster(file_path, tuples, out_path=out_file, compress=True)

                cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file])
                container_md = __internal__.merge_container_md(container_md, cur_md)

            elif filename.endswith('.las'):
                # If file is LAS, we can merge with any existing scan+plot output safely
                out_path = os.path.join(check_md['working_folder'], plot_name)
                out_file = os.path.join(out_path, filename)
                merged_out = os.path.join(out_path, os.path.splitext(os.path.basename(filename))[0] + '_merged.las')
                merged_txt = merged_out.replace('.las', '_contents.txt')
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

                if not __internal__.check_already_merged(merged_txt, file_path):
                    __internal__.clip_las(file_path, tuples, out_path=out_file, merged_path=merged_out)
                    with open(merged_txt, 'a') as contents:
                        contents.write(file_path + "\n")

                    cur_md = __internal__.prepare_container_md(plot_name, plot_md, sensor, file_path, [out_file, merged_out])
                    container_md = __internal__.merge_container_md(container_md, cur_md)
                else:
                    logging.info("Skipping already merged LAS data")

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
