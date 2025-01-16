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
        self.tools = [
            CreateLiDARProducts,
            CalculateOptimalPath,
            CreateLineProfile,
            CreateLiDARClassCostRaster,
            CalculateTopographicPositionalIndex,
            ]


class CreateLiDARProducts(object):
    def __init__(self):
        self.label = "Create LiDAR Products"
        self.description = "Create LiDAR Products"

    def getParameterInfo(self):

        # Input LAS file
        input_LAS_Files=arcpy.Parameter(
            displayName="Input LAS/LAZ File(s)",
            name="Input_LAS_Files",
            datatype="DEFile",
            multiValue=True,
            )
        
        # Convert LAZ
        extract_laz=arcpy.Parameter(
            displayName="Extract LAZ File(s)",
            name="Extract_LAZ",
            datatype="Boolean",
            parameterType="Optional",
            enabled=False,
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
        spatial_ref=arcpy.Parameter(
            displayName="Spatial Reference",
            name="Spatial_Reference",
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
        
        classify_height=arcpy.Parameter(
            displayName="Reclassify LAS by Height",
            name="Classify_LAS_Height",
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
        
        classify_overlap=arcpy.Parameter(
            displayName="Reclassify LAS Overlap",
            name="Classify_LAS_Overlap",
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

        # Sinks
        fill_sinks=arcpy.Parameter(
            displayName="Fill Sinks",
            name="Fill Sinks",
            datatype="Boolean",
            parameterType="Optional",
            category="Surfaces",
            enabled=False,
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
        # Classification
        out_classification=arcpy.Parameter(
            displayName="Classification",
            name="Create_Classification",
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
            enabled=False,
            )
        aggregate_polygons=arcpy.Parameter(
            displayName="Aggregate Polygons (remove holes and small features)",
            name="Agggregate_Polygons",
            datatype="Boolean",
            parameterType="Optional",
            category='Polygons',
            enabled=False,
            )
        minimum_area=arcpy.Parameter(
            displayName="Minimum Polygon/Hole Area (square meters)",
            name="Minimum_Area",
            datatype="Long",
            parameterType="Optional",
            category='Polygons',
            enabled=False,
            )

        return [
            input_LAS_Files, extract_laz, output_folder, out_name,
            spatial_ref, extent_polygon,
            classify_building, classify_ground, classify_height,
            classify_noise, classify_overlap,
            out_dem, fill_sinks, 
            out_slope, out_hillshade, out_intensity, out_classification,
            out_contours, out_contours_interval,
            polygons, aggregate_polygons, minimum_area]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [input_LAS_Files, extract_laz, output_folder, out_name,
         spatial_ref, extent_polygon,
         classify_building, classify_ground, classify_height,
         classify_noise, classify_overlap,
         out_dem, fill_sinks, 
         out_slope, out_hillshade, out_intensity, out_classification,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

        extent_polygon.filter.list = ["Polygon"]
        
        if input_LAS_Files.value and not spatial_ref.altered:
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            las_inputs = [i for i in inputs if os.path.splitext(i)[1].lower() == '.las']
            srs = [arcpy.Describe(i).spatialReference for i in las_inputs if i]
            if srs:
                spatial_ref.value = srs[0]

        if input_LAS_Files.value:             
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]   
            # Can't update the classification of points in zlas/laz files..
            if any([i for i in inputs if os.path.splitext(i)[1].lower() == '.laz']):
                extract_laz.enabled = True
                if not _3D or not extract_laz.value:
                    classify_building.value = False
                    classify_building.enabled = False
                    classify_ground.value = False
                    classify_ground.enabled = False
                    classify_height.enabled = False
                    classify_noise.value = False
                    classify_noise.enabled = False
                    classify_overlap.enabled = False
                else:
                    classify_building.enabled = True
                    classify_ground.enabled = True
                    classify_height.enabled = True
                    classify_noise.enabled = True
                    classify_overlap.enabled = True
            else:
                extract_laz.enabled = False
                classify_building.enabled = True
                classify_ground.enabled = True
                classify_height.enabled = True
                classify_noise.enabled = True
                classify_overlap.enabled = True
                    
        if out_dem.value:
            fill_sinks.enabled = True
            out_slope.enabled = True
            out_hillshade.enabled = True
            out_contours.enabled = True
        else:
            fill_sinks.value = None
            fill_sinks.enabled = False
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
            
        if out_classification.value:
            polygons.enabled = True
        else:
            polygons.enabled = False
            
        polygons.filter.type = "ValueList"
        polygons.filter.list = ['Low Vegetation', 'Medium Vegetation', 'High Vegetation', 'Building', 'Low Point',
                                'Water', 'Rail', 'Road Surface', 'Bridge Deck','Transmission Tower', 'Unassigned'
                                ]
        
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
        [input_LAS_Files, extract_laz, output_folder, out_name,
         spatial_ref, extent_polygon,
         classify_building, classify_ground, classify_height,
         classify_noise, classify_overlap,
         out_dem, fill_sinks, 
         out_slope, out_hillshade, out_intensity, out_classification,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

        # Set error on non-(.las/.lasz) input
        if input_LAS_Files.value:
            compressed_size = 0
            inputs = [i.strip("'") for i in input_LAS_Files.valueAsText.split(';')]
            for i in inputs:
                ext = os.path.splitext(i)[1].lower()
                if not ext in ALLOWED_FORMATS:
                    input_LAS_Files.setErrorMessage(f'Input: [{i}] is not a las/zlas/laz file..')
                elif ext == ".laz":
                    try:
                        compressed_size += os.stat(i).st_size
                    except: pass
                    
            # Set error on conflicting spatial references
            las_inputs = [i for i in inputs if os.path.splitext(i)[1].lower() == '.las']
            srs = [arcpy.Describe(i).spatialReference.name for i in las_inputs]
            if len(set(srs)) > 1:
                input_LAS_Files.setErrorMessage(f'Inputs must be in the same spatial reference..\n{srs}')

            # Set error on duplicate inputs
            if len(set(inputs)) != len(inputs):
                    err_mg = ("Duplicate input data: the same file cannnot be used more than once.")
                    input_LAS_Files.setErrorMessage(err_mg)
        
            if extract_laz.value:
                extract_laz.setWarningMessage(
                    f"Compressed Size: {round(compressed_size / 1073741824, 2)} gb"
                    )
        
        return

    def execute(self, parameters, messages):
        [input_LAS_Files, extract_laz, output_folder, out_name,
         spatial_ref, extent_polygon,
         classify_building, classify_ground, classify_height,
         classify_noise, classify_overlap,
         out_dem, fill_sinks, 
         out_slope, out_hillshade, out_intensity, out_classification,
         out_contours, out_contours_interval,
         polygons, aggregate_polygons, minimum_area] = parameters

        try:
            
            # Optionally extract laz
            if extract_laz.value:
                for laz_file in [f for f in input_LAS_Files.values if str(f).lower().endswith(".laz")]:
                    arcpy.conversion.ConvertLas(
                        in_las=laz_file, target_folder=os.path.dirname(str(laz_file))
                        )
                input_files = [str(f).lower().replace(".laz", ".las") for f in input_LAS_Files.values]
            else:
                input_files = input_LAS_Files.values
            
            # Create a geodatabase
            arcpy.CreateFileGDB_management(output_folder.value, out_name.value)

            # Set workspace to fGDB
            out_gdb = os.path.join(output_folder.valueAsText, out_name.valueAsText + '.gdb')
            arcpy.env.workspace = out_gdb

            # Create Las Dataset
            arcpy.AddMessage('Creating LAS Dataset..')
            out_lasd = os.path.join(output_folder.valueAsText, out_name.valueAsText + '.lasd')
            arcpy.management.CreateLasDataset(
                input=input_files,
                out_las_dataset=out_lasd,
                folder_recursion='NO_RECURSION',
                in_surface_constraints=None,
                spatial_reference=spatial_ref.value,
                compute_stats="COMPUTE_STATS",
                relative_paths="ABSOLUTE_PATHS",
                )
            
            # Update LAS classifications
            if classify_ground.value == 1:
                try:
                    arcpy.AddMessage('Updating Ground Classifications..')
                    arcpy.ddd.ClassifyLasGround(out_lasd, "STANDARD", "REUSE_GROUND", "0.5 Meters")
                except Exception as e:
                    arcpy.AddWarning(f'Could not update Ground Classifications: {e}')   
            if classify_building.value == 1:
                try:
                    arcpy.AddMessage('Updating Building classifications..')
                    arcpy.ddd.ClassifyLasBuilding(out_lasd, "2 Meters", "4 SquareMeters")
                except Exception as e:
                    arcpy.AddWarning(f'Could not update Building Classifications: {e}')   
            if classify_height.value == 1:
                try:
                    arcpy.AddMessage('Updating Height Classifications..')
                    arcpy.ddd.ClassifyLasByHeight(out_lasd, "GROUND", "3 5;4 25;5 50")
                except Exception as e:
                    arcpy.AddWarning(f'Could not update Height Classifications: {e}')   
            if classify_noise.value == 1:
                try:
                    arcpy.AddMessage('Updating Noise Classifications..')
                    arcpy.ddd.ClassifyLasNoise(out_lasd)
                except Exception as e:
                    arcpy.AddWarning(f'Could not update Noise Classifications: {e}')   
            if classify_overlap.value == 1:
                try:
                    arcpy.AddMessage('Updating Overlap Classifications..')
                    arcpy.ddd.ClassifyLasOverlap(out_lasd)
                except Exception as e:
                    arcpy.AddWarning(f'Could not update Overlap Classifications: {e}')   
            
           
            # Create DEM and other surfaces
            if out_dem.value:
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
                    overlap="EXCLUDE_OVERLAP",
                    )
                arcpy.AddMessage('Creating DEM Surface..')
                dem = arcpy.conversion.LasDatasetToRaster(
                    dem_lasd, os.path.join(out_gdb, "Elevation"), "ELEVATION", "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", 1, 1
                )
                if extent_polygon.value:
                    outExtractByMask = arcpy.sa.ExtractByMask(dem, extent_polygon.value)
                    outExtractByMask.save(os.path.join(out_gdb, "Elevation"))
                arcpy.management.Delete(dem_lasd)
            
            # Fill sinks
            if fill_sinks.value:
                try:
                    arcpy.AddMessage('Filling Sinks in DEM..')
                    filled_dem = arcpy.sa.Fill(dem)
                    filled_dem.save(os.path.join(out_gdb, "Elevation"))
                except Exception as e:
                    arcpy.AddWarning(f'Could not fill sinks in DEM: {e}')
                    
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

            if out_classification.value:
                try:
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
    
                    # Create attribute table
                    arcpy.management.BuildRasterAttributeTable(out_raster)

                except Exception as e:
                    arcpy.AddWarning(f'Could not create classified surface: {e}')
                     
            if out_contours.value:
                try:
                    arcpy.AddMessage('Creating Contours..')
                    try:
                        units = arcpy.Describe(dem).spatialReference.VCS.linearUnitName
                        contour_name = f"Contours_{out_contours_interval.value}_{units.replace(' ', '_')}"
                    except:
                        contour_name = f"Contours_{out_contours_interval.value}"
                    contours = os.path.join(out_gdb, contour_name)
                    arcpy.sa.Contour(dem, contours, out_contours_interval.value)
                except Exception as e:
                    arcpy.AddWarning(f'Could not create contours: {e}')    
                        
            # Extract Classes as polygons
            if out_classification.value and polygons.value and arcpy.Exists(class_raster):
                try:
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

                except Exception as e:
                    arcpy.AddWarning(f'Could not create polygons: {e}')
                    
                arcpy.AddMessage('Done..')

        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return


class CalculateOptimalPath(object):
    def __init__(self):
        self.label = "Calculate Optimal Path"
        self.description = "Calculate Optimal Path"
        self.category = "Analysis"

    def getParameterInfo(self):

        # Source
        source=arcpy.Parameter(
            displayName="Source Feature",
            name="Source",
            datatype="GPFeatureLayer",
            )

        # Destination
        destination=arcpy.Parameter(
            displayName="Destination Feature",
            name="Destination",
            datatype="GPFeatureLayer",
            )

        # DEM
        elevation=arcpy.Parameter(
            displayName="Elevation Raster",
            name="Elevation",
            datatype=["DERasterDataset", "GPRasterLayer"],
            )

        # Costs
        costs=arcpy.Parameter(
            displayName="Cost Raster",
            name="Costs",
            datatype=["DERasterDataset", "GPRasterLayer"],
            )
                
        # Barriers
        barriers=arcpy.Parameter(
            displayName="Barrier Features",
            name="Barriers",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            )

        # Output
        optimal_path=arcpy.Parameter(
            displayName="Output Optimal Path",
            name="Output Optimal Path",
            datatype="DEFeatureClass",
            direction="Output",
            )

        return [source, destination, elevation, costs, barriers, optimal_path]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [source, destination, elevation, costs, barriers, optimal_path] = parameters
            
        return

    def updateMessages(self, parameters):
        [source, destination, elevation, costs, barriers, optimal_path] = parameters     

        return

    def execute(self, parameters, messages):
        [source, destination, elevation, costs, barriers, optimal_path] = parameters

        try:
            
            arcpy.AddMessage("Creating distance accumulation and back direction rasters..")
            out_distance_allocation_raster = arcpy.sa.DistanceAllocation(
                in_source_data=source.value,
                in_barrier_data=barriers.value,
                in_surface_raster=elevation.value,
                in_cost_raster=costs.value,
                out_distance_accumulation_raster="distance_accumulation",
                out_back_direction_raster="back_direction",
                )
            out_distance_allocation_raster.save("distance_allocation")

            arcpy.AddMessage("Calculating optimal path..")
            arcpy.sa.OptimalPathAsLine(
                in_destination_data=destination.value,
                in_distance_accumulation_raster="distance_accumulation",
                in_back_direction_raster="back_direction",
                out_polyline_features=optimal_path.value,
                path_type="EACH_ZONE",
            )

            for item in ["distance_accumulation", "back_direction", "distance_allocation"]:
                try:
                    arcpy.management.Delete(item)
                except: pass

        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
    

class CreateLineProfile(object):
    def __init__(self):
        self.label = "Create 3D Line Profile"
        self.description = "Create 3D Line Profile"
        self.category = "Analysis"

    def getParameterInfo(self):

        # Source Line
        source=arcpy.Parameter(
            displayName="Source Line Feature",
            name="Source_Line",
            datatype="Feature Layer",
            )

        # DEM
        elevation=arcpy.Parameter(
            displayName="Elevation Raster",
            name="Elevation",
            datatype=["DERasterDataset", "GPRasterLayer"],
            )

        # Method
        method=arcpy.Parameter(
            displayName="Sampling Method",
            name="Method",
            datatype="GPString"
            )
                
        # Interval
        interval=arcpy.Parameter(
            displayName="Distance Interval",
            name="Interval",
            datatype="GPLong",
            parameterType="Optional",
            enabled=False,
            )

        # Units
        units=arcpy.Parameter(
            displayName="Distance Unit",
            name="Units",
            datatype="GPString",
            parameterType="Optional",
            enabled=False,
            )

        # Percent
        percent=arcpy.Parameter(
            displayName="Sampling Percent",
            name="Percent",
            datatype="GPLong",
            parameterType="Optional",
            enabled=False,
            )        

        # Output
        output_line=arcpy.Parameter(
            displayName="Output 3D Line Profile",
            name="Output 3D Line Profile",
            datatype="DEFeatureClass",
            direction="Output",
            )

        return [source, elevation, method, interval, units, percent, output_line]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [source, elevation, method, interval, units, percent, output_line] = parameters
        
        # Enfore line
        source.filter.list = ["Polyline"]
                
        # Filter Method Options
        method.filter.type = "ValueList"
        method.filter.list = ["Distance", "Percentage"]
        
        units.filter.type = "ValueList"
        units.filter.list = ["Feet", "Meters"]
        
        # Filter method paramters
        if method.value:
            if method.value == "Distance":
                interval.enabled = True
                if not interval.altered: interval.value = 10
                units.enabled = True
                if not units.altered: units.value = "Feet"
                percent.value = None; percent.enabled = False
            else: #Percent
                percent.enabled = True
                if not percent.altered: percent.value = 5
                interval.value = None; interval.enabled = False
                units.value = None; units.enabled = False            
                
        return

    def updateMessages(self, parameters):
        [source, elevation, method, interval, units, percent, output_line] = parameters
  
        # Enforce Single-part
        if source.value:
            cnts = int(arcpy.management.GetCount(source.value).getOutput(0))
            if cnts != 1:
                sels = [i for i in arcpy.Describe(source.value).FIDSet.split(";") if i]
                if len(sels) != 1:
                    source.setErrorMessage('Select a single feature from the input layer..')
        return

    def execute(self, parameters, messages):
        [source, elevation, method, interval, units, percent, output_line] = parameters

        try:
            # Create points along line
            arcpy.AddMessage('Generating sample points..')
            pts = arcpy.management.GeneratePointsAlongLines(
                Input_Features=source.value,
                Output_Feature_Class=r"points",
                Point_Placement=method.value.upper(),
                Distance=f"{interval.value} {units.value}" if method.valueAsText == "Distance" else None,
                Percentage=percent.value if method.valueAsText == "Percentage" else 0,
                Include_End_Points="END_POINTS",
                Add_Chainage_Fields="ADD_CHAINAGE",
            )
            # Extract z values to points
            arcpy.AddMessage('Extracting elevation values..')
            ptsz = arcpy.sa.ExtractValuesToPoints(
                in_point_features=pts,
                in_raster=elevation.value,
                out_point_features=r"points_z",
                interpolate_values="INTERPOLATE",
                add_attributes="VALUE_ONLY"
            )
            # Convert points to 3D
            arcpy.AddMessage('Converting sample points to 3D..')
            pts3d = arcpy.ddd.FeatureTo3DByAttribute(
                in_features=ptsz,
                out_feature_class=r"points_3d",
                height_field="RASTERVALU",
            )
            # Convert 3D points back to a line
            arcpy.AddMessage('Generating 3D line..')
            arcpy.management.PointsToLine(
                Input_Features=pts3d,
                Output_Feature_Class=output_line.value,
                Line_Field=None,
                Sort_Field="ORIG_SEQ",
                Close_Line="NO_CLOSE",
                Line_Construction_Method="CONTINUOUS",
                Attribute_Source="NONE",
                Transfer_Fields=None
            )
        
            for item in [pts, ptsz, pts3d]:
                arcpy.management.Delete(item)
                
        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
    
    
class CreateLiDARClassCostRaster(object):
    def __init__(self):
        self.label = "Create LiDAR Class Cost Raster"
        self.description = "Create LiDAR Class Cost Raster"
        self.category = "Analysis"

    def getParameterInfo(self):

        # Source Raster
        source=arcpy.Parameter(
            displayName="Source Raster",
            name="Source_Raster",
            datatype="DERasterDataset",
            )
        # Class Costs
        costs=arcpy.Parameter(
            displayName="Class Costs",
            name="Class Costs",
            datatype="GPValueTable",
            )
        costs.columns = [['String', 'Class'], ['Long', 'Cost']]

        return [source, costs]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [source, costs] = parameters           
        
        if not costs.altered:
            costs.values = [
                [k, 1] for k in LiDAR_CLASS_CODES.values()
            ]
            
        return

    def updateMessages(self, parameters):
        [source, costs] = parameters
        return

    def execute(self, parameters, messages):
        [source, costs] = parameters

        try:
            class_cost_dict = {class_name: cost for class_name, cost in costs.values}
            arcpy.management.BuildRasterAttributeTable(source.value)
            for field_name, field_type in (("Name", "Text"), ("Cost", "Long")):
                try:
                    arcpy.management.AddField(source.value, field_name=field_name, field_type=field_type)
                except: pass
            with arcpy.da.UpdateCursor(source.value, ["Value", "Name", "Cost"]) as cur:
                for value, name, cost in cur:
                    if value:
                        name = LiDAR_CLASS_CODES.get(value, "Other")
                        cost = class_cost_dict.get(name, 1)
                        cur.updateRow((value, name, cost))
                    
        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
    
    
class CalculateTopographicPositionalIndex(object):
    def __init__(self):
        self.label = "Calculate Topographic Positional Index"
        self.description = "Calculate Topographic Positional Index"
        self.category = "Analysis"

    def getParameterInfo(self):

        # DEM
        elevation=arcpy.Parameter(
            displayName="Elevation Raster",
            name="Elevation",
            datatype=["DERasterDataset", "GPRasterLayer"],
            )
        shape=arcpy.Parameter(
            displayName="Neighborhood Shape",
            name="Neighborhood_Shape",
            datatype="GPString",
            )
        distance=arcpy.Parameter(
            displayName="Neighborhood Distance",
            name="Neighborhood_Distance",
            datatype="GPLong",
        )
        units=arcpy.Parameter(
            displayName="Neighborhood Units",
            name="Neighborhood_Units",
            datatype="GPString",
            )
        output=arcpy.Parameter(
            displayName="Output Raster",
            name="Output_Raster",
            datatype="DERasterDataset",
            direction="Output",
            )
        
        return [elevation, shape, distance, units, output]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        [elevation, shape, distance, units, output] = parameters
        
        shape.filter.type = "ValueList"
        shape.filter.list = ["Circle", "Rectangle"]
        if not shape.altered: shape.value = "Rectangle"
        units.filter.type = "ValueList"
        units.filter.list = ["Cell", "Map"]
        if not units.altered: units.value = "Cell"
        if not distance.altered: distance.value = 10
        return

    def updateMessages(self, parameters):
        [elevation, shape, distance, units, output] = parameters
        return

    def execute(self, parameters, messages):
        [elevation, shape, distance, units, output] = parameters
        
        try:
            shapes = {
                    "Circle": f"{shape.value} {distance.value} {units.value}",
                    "Rectangle": f"{shape.value} {distance.value} {distance.value} {units.value}",
                } 
            mean = arcpy.ia.FocalStatistics(elevation.valueAsText, shapes[shape.valueAsText], "MEAN")
            diff = arcpy.ddd.Minus(elevation.value, mean, output.valueAsText)
  
        except:
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        return
