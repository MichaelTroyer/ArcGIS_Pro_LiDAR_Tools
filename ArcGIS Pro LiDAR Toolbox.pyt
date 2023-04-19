# -*- coding: utf-8 -*-

"""
Michael Troyer

michael.troyer@usda.gov


Summary:
LiDAR streamlining


"""
import os
import traceback
import arcpy

arcpy.env.overwriteOutput = True


LiDAR_CLASS_CODES = {
    0  : 'Never classified',
    1  : 'Unassigned',
    2  : 'Ground',
    3  : 'Low Vegetation',
    4  : 'Medium Vegetation',
    5  : 'High Vegetation',
    6  : 'Building',
    7  : 'Low Point',
    8  : 'Reserved',
    9  : 'Water',
    10 : 'Rail',
    11 : 'Road Surface',
    12 : 'Reserved',
    13 : 'Wire - Guard (Shield)',
    14 : 'Wire - Conductor (Phase)',
    15 : 'Transmission Tower',
    16 : 'Wire-Structure Connector (Insulator)',
    17 : 'Bridge Deck',
    18 : 'High Noise',
    }
REVERSE_LOOKUP = {v:k for k, v in LiDAR_CLASS_CODES.items()}


def build_where_clause(table, field, valueList):
    """Takes a list of values and constructs a SQL WHERE
    clause to select those values within a given field and table."""
    fieldDelimited = arcpy.AddFieldDelimiters(arcpy.Describe(table).path, field)
    fieldType = arcpy.ListFields(table, field)[0].type
    if str(fieldType) == 'String':
        valueList = ["'%s'" % value for value in valueList]
    whereClause = "%s IN(%s)" % (fieldDelimited, ', '.join(map(str, valueList)))
    return whereClause


class Toolbox(object):
    def __init__(self):
        self.label = "LiDAR Tools"
        self.alias = "LiDAR_Tools"
        self.tools = [CreateLiDARProducts]


class CreateLiDARProducts(object):
    def __init__(self):
        self.label = "Create LiDAR Products"
        self.description = "Create LiDAR Products"

    def getParameterInfo(self):

        # Input LAZ file
        input_LAS_Files=arcpy.Parameter(
            displayName="Input LAS/ZLAS File(s)",
            name="Input_LAS_Files",
            datatype="DEFile",
            multiValue=True,
            )

        # Output location
        output_folder=arcpy.Parameter(
            displayName="Output Folder",
            name="Output_Folder",
            datatype="DEFolder",
            )

        # Output name
        out_name = arcpy.Parameter(
            displayName="Output Name",
            name="Output_Name",
            datatype="GPString",
            )

        # Projection
        projection=arcpy.Parameter(
            displayName="LAS File Projection",
            name="LAS_File_Projection",
            datatype="GPSpatialReference",
            )
                
        # Extent
        extent_polygon=arcpy.Parameter(
            displayName="Extent Polygon",
            name="Extent_Polygon",
            datatype="Feature Layer",
            parameterType="Optional",
            )

        # DEM
        out_dem=arcpy.Parameter(
            displayName="Digital Elevation Model",
            name="create_dem",
            datatype="Boolean",
            parameterType="Optional",
            category='Rasters',
            )
        
        # Slope
        out_slope=arcpy.Parameter(
            displayName="Slope",
            name="Create_Slope",
            datatype="Boolean",
            parameterType="Optional",
            category='Rasters',
            enabled=False,
            )
        
        # Hillshade
        out_hillshade=arcpy.Parameter(
            displayName="Hillshade",
            name="Create_Hillshade",
            datatype="Boolean",
            parameterType="Optional",
            category='Rasters',
            enabled=False,
            )

        # Intensity
        out_intensity=arcpy.Parameter(
            displayName="Return Intensity",
            name="Create_Intensity",
            datatype="Boolean",
            parameterType="Optional",
            category='Rasters',
            )

        # Classifications
        classify_building=arcpy.Parameter(
            displayName="Reclassify LAS Buildings",
            name="Classify_LAS_Buildings",
            datatype="Boolean",
            parameterType="Optional",
            category='Classification',
            )
        
        classify_ground=arcpy.Parameter(
            displayName="Reclassify LAS Ground",
            name="Classify_LAS_Ground",
            datatype="Boolean",
            parameterType="Optional",
            category='Classification',
            )
        
        classify_noise=arcpy.Parameter(
            displayName="Reclassify LAS Noise",
            name="Classify_LAS_Noise",
            datatype="Boolean",
            parameterType="Optional",
            category='Classification',
            )

        #Polygons
        polygons=arcpy.Parameter(
            displayName="Extract LiDAR Classes as Polygons",
            name="Extract_LiDAR_Polygons",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            category='Polygons',
            )

        return [
            input_LAS_Files, output_folder, out_name, projection, extent_polygon,
            classify_building, classify_ground, classify_noise,
            out_dem, out_slope, out_hillshade, out_intensity, 
            polygons]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [input_LAS_Files, output_folder, out_name, projection, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity, 
         polygons] = parameters

        extent_polygon.filter.list = ["Polygon"]
        if input_LAS_Files.value:
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            inputs = [i for i in inputs if os.path.splitext(i)[1].lower() in ('.las', '.zlas')] 
            srs = [arcpy.Describe(i).spatialReference for i in inputs]
            projection.value = srs[0]
            
            # Can't update the classification of points in zlas files..
            if any([i for i in inputs if os.path.splitext(i)[1].lower() == '.zlas']):
                classify_building.value = None
                classify_building.enabled = False
                classify_ground.value = None
                classify_ground.enabled = False
                classify_noise.value = None
                classify_noise.enabled = False
            else:
                classify_building.enabled = True
                classify_ground.enabled = True
                classify_noise.enabled = True
                
        if out_dem.value:
            out_slope.enabled = True
            out_hillshade.enabled = True
        else:
            out_slope.value = None
            out_slope.enabled = False
            out_hillshade.value = False
            out_hillshade.enabled = False

        polygons.filter.type = "ValueList"
        polygons.filter.list = ['Ground', 'Low Vegetation', 'Medium Vegetation', 'High Vegetation', 'Building',
                                'Low Point', 'Water', 'Rail', 'Road Surface', 'Bridge Deck',]
        return

    def updateMessages(self, parameters):
        [input_LAS_Files, output_folder, out_name, projection, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity, 
         polygons] = parameters

        # Set error on non-(.las/.lasz) input
        if input_LAS_Files.value:
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            for i in inputs:
                if not os.path.splitext(i)[1].lower() in ('.las', '.zlas'):
                    input_LAS_Files.setErrorMessage(f'Input: [{i}] is not a las/zlas file..')
                    
            # Set error on conflicting spatial references
            inputs = [i for i in inputs if os.path.splitext(i)[1].lower() in ('.las', '.zlas')]
            srs = [arcpy.Describe(i).spatialReference.name for i in inputs]
            if len(set(srs)) > 1:
                input_LAS_Files.setErrorMessage(f'Inputs must be in the same projection..\n{srs}')

        # # Verify a single polygon is selected
        # if extent_polygon.value:
        #     if not len([i for i in arcpy.Describe(extent_polygon.value).FIDSet.split(';')]) == 1:
        #         extent_polygon.setErrorMessage('Select a single polygon from the input layer..')
                
        return

    def execute(self, parameters, messages):
        [input_LAS_Files, output_folder, out_name, projection, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity, 
         polygons] = parameters

        try:
            # Create a geodatabase
            arcpy.CreateFileGDB_management(output_folder.value, out_name.value)

            # Set workspace to fGDB
            out_gdb = os.path.join(output_folder.valueAsText, out_name.valueAsText + '.gdb')
            arcpy.env.workspace = out_gdb

            # Create Las Dataset
            arcpy.AddMessage('Creating LAS Dataset..')
            out_lasd = os.path.join(output_folder.valueAsText, out_name.valueAsText + '.lasd')
            arcpy.management.CreateLasDataset(
                input=input_LAS_Files.values,
                out_las_dataset=out_lasd,
                folder_recursion='NO_RECURSION',
                in_surface_constraints=None,
                spatial_reference=projection.value,
                compute_stats="COMPUTE_STATS",
                relative_paths="ABSOLUTE_PATHS",
                create_las_prj="ALL_FILES",
                )
            
            # Update LAS classifications
            if classify_ground.value == 1:
                arcpy.AddMessage('Updating Ground Classifications..')
                arcpy.ddd.ClassifyLasGround(out_lasd, "STANDARD", "REUSE_GROUND", "0.5 Meters")

            if classify_building.value == 1:
                arcpy.AddMessage('Updating Building Classifications..')
                arcpy.ddd.ClassifyLasBuilding(out_lasd, "2 Meters", "4 SquareMeters")

            if classify_noise.value == 1:
                arcpy.AddMessage('Updating Noise Classifications..')
                arcpy.ddd.ClassifyLasNoise(
                    out_lasd, "ISOLATION", "CLASSIFY", "NO_WITHHELD", "COMPUTE_STATS", None, None, None, 10,
                    "8 Meters", "8 Meters", "DEFAULT", "PROCESS_EXTENT", None, "UPDATE_PYRAMID"
                    )

            # Make a bare-ground only lasd layer for DEM and subsequent calculations
            dem_lasd = arcpy.management.MakeLasDatasetLayer(
                in_las_dataset=out_lasd,
                out_layer='LasDatasetLayer',
                class_code="2",
                return_values=None,
                no_flag="INCLUDE_UNFLAGGED",
                synthetic="INCLUDE_SYNTHETIC",
                keypoint="INCLUDE_KEYPOINT",
                withheld="EXCLUDE_WITHHELD",
                surface_constraints=None,
                overlap="INCLUDE_OVERLAP",
                )
            
            # Create DEM and Intensity surfaces
            if out_dem.value == 1:
                arcpy.AddMessage('Creating DEM Surface..')
                dem = arcpy.conversion.LasDatasetToRaster(
                    dem_lasd, os.path.join(out_gdb, "Elevation"), "ELEVATION", "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", 1, 1
                )
                if extent_polygon.value:
                    outExtractByMask = arcpy.sa.ExtractByMask(dem, extent_polygon.value)
                    outExtractByMask.save(os.path.join(out_gdb, "Elevation"))
            
            if out_slope.value:
                arcpy.AddMessage('Creating Slope Surface..')
                slope = arcpy.ddd.Slope(dem, os.path.join(out_gdb, "Slope"), "DEGREE", 1, "PLANAR", "METER")
                    
            if out_hillshade.value:
                arcpy.AddMessage('Creating Hillshade..')
                shade = arcpy.ddd.HillShade(dem, os.path.join(out_gdb, "Hillshade"), 315, 45, "NO_SHADOWS", 1)
                    
            if out_intensity.value == 1:
                arcpy.AddMessage('Creating Intensity Surface..')
                intensity = arcpy.conversion.LasDatasetToRaster(
                    out_lasd, os.path.join(out_gdb, "Intensity"), "INTENSITY", "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", 1, 1
                )
                if extent_polygon.value:
                    outExtractByMask = arcpy.sa.ExtractByMask(intensity, extent_polygon.value)
                    outExtractByMask.save(os.path.join(out_gdb, "Intensity"))

            # Extract Classes as polygons
            if polygons.value:
                arcpy.AddMessage('Generating LAS Class Statistics..')
                class_raster = os.path.join(out_gdb, "Classes")
                out_raster = arcpy.management.LasPointStatsAsRaster(
                    in_las_dataset=out_lasd,
                    out_raster=class_raster,
                    method='PREDOMINANT_CLASS',
                    sampling_type='CELLSIZE',
                    sampling_value=1,
                    )
                if extent_polygon.value:
                    outExtractByMask = arcpy.sa.ExtractByMask(class_raster, extent_polygon.value)
                    outExtractByMask.save(os.path.join(out_gdb, "Classes"))
                    
                # Raster to polygon
                arcpy.AddMessage('Creating Class Polygons..')
                class_polygons = os.path.join(out_gdb, "Class_Polygons")
                arcpy.conversion.RasterToPolygon(
                    in_raster=class_raster,
                    out_polygon_features=class_polygons,
                    simplify='SIMPLIFY',
                    raster_field='VALUE',
                    create_multipart_features="SINGLE_OUTER_PART",
                    max_vertices_per_feature=None,
                    )

                arcpy.AddField_management(class_polygons, 'Class', 'String')

                # Select classes
                arcpy.AddMessage('Extracting Classes..')
                where = build_where_clause(class_polygons, 'gridcode', [REVERSE_LOOKUP[v.strip("'")] for v in polygons.valueAsText.split(';')])
                polygon_layer = arcpy.MakeFeatureLayer_management(class_polygons, 'Class_Polygons_lyr', where)

                with arcpy.da.UpdateCursor(polygon_layer, ['gridcode', 'Class']) as cur:
                    for row in cur:
                        try:
                            row[1] = LiDAR_CLASS_CODES[row[0]]
                            cur.updateRow(row)
                        except:
                            row[1] = 'Undefined'

                arcpy.CopyFeatures_management(polygon_layer, os.path.join(out_gdb, "Classified_Polygons"))

                arcpy.Delete_management(class_polygons)
                arcpy.Delete_management(polygon_layer)

                arcpy.AddMessage('Done..')

        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
