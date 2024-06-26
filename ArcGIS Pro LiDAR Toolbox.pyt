# -*- coding: utf-8 -*-

"""
Michael Troyer

michael.troyer@usda.gov


Summary:
LiDAR streamlining


"""
import math
import os
import traceback
import arcpy

arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("Spatial")

# Check for 3D license. If no 3D license, LiDAR point reclassification is not available. 
if arcpy.CheckExtension("3D") == "Available":
    arcpy.CheckOutExtension("3D")
    _3D = True
else:
    _3D = False


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

ALLOWED_FORMATS = ('.las', '.zlas', '.laz')

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
        self.alias = "LiDARTools"
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

        # Coordiante System
        coord_sys=arcpy.Parameter(
            displayName="LAS File Coordinate System",
            name="LAS_File_Coordinate_System",
            datatype="GPSpatialReference",
            )
                
        # Extent
        extent_polygon=arcpy.Parameter(
            displayName="Extent Polygon",
            name="Extent_Polygon",
            datatype="Feature Layer",
            parameterType="Optional",
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
        
        # DEM
        out_dem=arcpy.Parameter(
            displayName="Digital Elevation Model",
            name="create_dem",
            datatype="Boolean",
            parameterType="Optional",
            category="Surfaces",
            )
        
        # Slope
        out_slope=arcpy.Parameter(
            displayName="Slope",
            name="Create_Slope",
            datatype="Boolean",
            parameterType="Optional",
            category="Surfaces",
            enabled=False,
            )
        
        # Hillshade
        out_hillshade=arcpy.Parameter(
            displayName="Hillshade",
            name="Create_Hillshade",
            datatype="Boolean",
            parameterType="Optional",
            category="Surfaces",
            enabled=False,
            )

        # Intensity
        out_intensity=arcpy.Parameter(
            displayName="Intensity",
            name="Create_Intensity",
            datatype="Boolean",
            parameterType="Optional",
            category="Surfaces",
            )
        
        #Contours
        out_contours=arcpy.Parameter(
            displayName="Create Contours",
            name="Create_Contours",
            datatype="Boolean",
            parameterType="Optional",
            category="Contours",
            enabled=False,
            )
        out_contours_interval=arcpy.Parameter(
            displayName="Contour Interval (in source units)",
            name="Contour_Interval",
            datatype="Long",
            parameterType="Optional",
            category="Contours",
            enabled=False,
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
        aggregate_polygons=arcpy.Parameter(
            displayName="Aggregate Polygons (remove holes and small features)",
            name="Agggregate_Polygons",
            datatype="Boolean",
            parameterType="Optional",
            category='Polygons',
            )
        minimum_area=arcpy.Parameter(
            displayName="Minimum Polygon/Hole Area (square meters)",
            name="Minimum_Area",
            datatype="Long",
            parameterType="Optional",
            category='Polygons',
            )

        return [
            input_LAS_Files, output_folder, out_name,
            coord_sys, extent_polygon,
            classify_building, classify_ground, classify_noise,
            out_dem, out_slope, out_hillshade, out_intensity,
            out_contours, out_contours_interval,
            polygons, aggregate_polygons, minimum_area]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [input_LAS_Files, output_folder, out_name,
         coord_sys, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

        extent_polygon.filter.list = ["Polygon"]
        
        if input_LAS_Files.value and not coord_sys.altered:
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            las_inputs = [i for i in inputs if os.path.splitext(i)[1].lower() == '.las']
            srs = [arcpy.Describe(i).spatialReference for i in las_inputs]
            coord_sys.value = srs[0]
        elif input_LAS_Files.value:             
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]   
            # Can't update the classification of points in zlas/laz files..
            if any([i for i in inputs if os.path.splitext(i)[1].lower() in ('.zlas', '.laz')]) or not _3D:
                classify_building.value = False
                classify_building.enabled = False
                classify_ground.value = False
                classify_ground.enabled = False
                classify_noise.value = False
                classify_noise.enabled = False
            else:
                classify_building.enabled = True
                classify_ground.enabled = True
                classify_noise.enabled = True
                
        if out_dem.value:
            out_slope.enabled = True
            out_hillshade.enabled = True
            out_contours.enabled = True
        else:
            out_slope.value = False
            out_slope.enabled = False
            out_hillshade.value = False
            out_hillshade.enabled = False
            out_contours.value = False
            out_contours.enabled = False
            
        if out_contours.value:
            out_contours_interval.enabled = True
            if not out_contours_interval.altered:
                out_contours_interval.value = 3
        else:
            out_contours_interval.value = None
            out_contours_interval.enabled = False

        polygons.filter.type = "ValueList"
        polygons.filter.list = ['Ground', 'Low Vegetation', 'Medium Vegetation', 'High Vegetation', 'Building',
                                'Low Point', 'Water', 'Rail', 'Road Surface', 'Bridge Deck','Transmission Tower',
                                'High Noise', 'Unassigned']
        
        if polygons.value:
            aggregate_polygons.enabled = True
        else:
            aggregate_polygons.value = None
            minimum_area.value = None
            aggregate_polygons.enabled = False
            minimum_area.enabled = False
            
        if aggregate_polygons.value:
            minimum_area.enabled = True
            if not minimum_area.value:
                minimum_area.value = 5
        else:
            minimum_area.value = None
            minimum_area.enabled = False
            
        return

    def updateMessages(self, parameters):
        [input_LAS_Files, output_folder, out_name,
         coord_sys, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

        # Set error on non-(.las/.lasz) input
        if input_LAS_Files.value:
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            for i in inputs:
                if not os.path.splitext(i)[1].lower() in ALLOWED_FORMATS:
                    input_LAS_Files.setErrorMessage(f'Input: [{i}] is not a las/zlas/laz file..')
                    
            # Set error on conflicting spatial references
            las_inputs = [i for i in inputs if os.path.splitext(i)[1].lower() == '.las']
            srs = [arcpy.Describe(i).spatialReference.name for i in las_inputs]
            if len(set(srs)) > 1:
                input_LAS_Files.setErrorMessage(f'Inputs must be in the same coordinate system..\n{srs}')

            # Set error on duplicate inputs
            if len(set(inputs)) != len(inputs):
                    err_mg = ("Duplicate input data: the same file cannnot be used more than once.")
                    input_LAS_Files.setErrorMessage(err_mg)
        
        
        return

    def execute(self, parameters, messages):
        [input_LAS_Files, output_folder, out_name,
         coord_sys, extent_polygon,
         classify_building, classify_ground, classify_noise,
         out_dem, out_slope, out_hillshade, out_intensity,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

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
                spatial_reference=coord_sys.value,
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
                ###
                # overlap="INCLUDE_OVERLAP",
                overlap="EXCLUDE_OVERLAP",
                ###
                )
            
            # Create DEM and Intensity surfaces
            if out_dem.value:
                arcpy.AddMessage('Creating DEM Surface..')
                dem = arcpy.conversion.LasDatasetToRaster(
                    dem_lasd, os.path.join(out_gdb, "Elevation"), "ELEVATION", "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", 1, 1
                )
                if extent_polygon.value:
                    outExtractByMask = arcpy.sa.ExtractByMask(dem, extent_polygon.value)
                    outExtractByMask.save(os.path.join(out_gdb, "Elevation"))
            
            if out_slope.value:
                try:
                    arcpy.AddMessage('Creating Slope Surface..')
                    slope = arcpy.sa.SurfaceParameters(dem, 'SLOPE', output_slope_measurement='Degree')
                    slope.save(os.path.join(out_gdb, "Slope"))
                except Exception as e:
                    arcpy.AddWarning(f'Could not create slope raster: {e}')
                
            if out_hillshade.value:
                try:
                    arcpy.AddMessage('Creating Hillshade..')
                    shade = arcpy.sa.Hillshade(dem); shade.save(os.path.join(out_gdb, "Hillshade"))
                except Exception as e:
                    arcpy.AddWarning(f'Could not create hillshade raster: {e}')                   

            if out_intensity.value:
                try:
                    arcpy.AddMessage('Creating Intensity Surface..')
                    intensity = arcpy.conversion.LasDatasetToRaster(
                        out_lasd, os.path.join(out_gdb, "Intensity"), "INTENSITY", "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", 1, 1
                    )
                    if extent_polygon.value:
                        outExtractByMask = arcpy.sa.ExtractByMask(intensity, extent_polygon.value)
                        outExtractByMask.save(os.path.join(out_gdb, "Intensity"))
                except Exception as e:
                    arcpy.AddWarning(f'Could not create intensity surface: {e}')

            if out_contours.value:
                try:
                    arcpy.AddMessage('Creating Contours..')
                    # TODO: allow alternative output units
                    units = arcpy.Describe(dem).spatialReference.VCS.linearUnitName
                    contours = os.path.join(
                        out_gdb, f"Contours_{out_contours_interval.value}_{units.replace(' ', '_')}")
                    arcpy.sa.Contour(dem, contours, out_contours_interval.value)
                except Exception as e:
                    arcpy.AddWarning(f'Could not create contours: {e}')    
            
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

                if aggregate_polygons.value:
                    arcpy.AddMessage('Aggregating Polygons..')
                    arcpy.cartography.AggregatePolygons(
                        in_features=polygon_layer,
                        out_feature_class=os.path.join(out_gdb, "Classified_Polygons"),
                        aggregation_distance=f"{round(math.sqrt(minimum_area.value))} Meters",
                        minimum_area=f"{minimum_area.value} SquareMeters",
                        minimum_hole_size=f"{minimum_area.value} SquareMeters",
                        aggregate_field="Class",
                        )
                else:
                    arcpy.CopyFeatures_management(polygon_layer, os.path.join(out_gdb, "Classified_Polygons"))

                arcpy.Delete_management(class_polygons)
                arcpy.Delete_management(polygon_layer)

                arcpy.AddMessage('Done..')

        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
