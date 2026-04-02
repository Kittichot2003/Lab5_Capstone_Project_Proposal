# Lab5_Capstone_Project_Proposal
# ==============================================================================
# PROJECT: AVATAR NA'VI HABITAT SUITABILITY MODEL (THAILAND) - LATEST VERSION
# ==============================================================================

# ---------------------------------------------------------
# 1. ติดตั้งไลบรารีและเชื่อมต่อระบบ
# ---------------------------------------------------------
import ee
import geemap

ee.Initialize(project='v6-dk-459909') 

# ---------------------------------------------------------
# 2. กำหนดขอบเขตและดึงข้อมูลพื้นฐาน (Base Data)
# ---------------------------------------------------------
Map = geemap.Map()
thailand = ee.FeatureCollection("FAO/GAUL/2015/level0").filter(ee.Filter.eq('ADM0_NAME', 'Thailand'))
Map.centerObject(thailand, 6)

start_date = '2020-01-01'
end_date = '2023-12-31'

# ความสูง (Elevation) และ ความลาดชัน (Slope)
elevation = ee.Image('USGS/SRTMGL1_003').clip(thailand)
elev_norm = elevation.unitScale(0, 2500).clamp(0, 1)

slope = ee.Terrain.slope(elevation)
slope_norm = slope.unitScale(0, 30).clamp(0, 1)

# สภาพแวดล้อม (NDVI, Rainfall, Temperature)
ndvi = ee.ImageCollection('MODIS/061/MOD13A1').filterDate(start_date, end_date).select('NDVI').mean().multiply(0.0001).clip(thailand)
ndvi_norm = ndvi.unitScale(0, 0.9).clamp(0, 1)

rainfall = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterDate(start_date, end_date).sum().divide(4).clip(thailand)
rain_norm = rainfall.unitScale(1000, 4000).clamp(0, 1)

temp = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR").filterDate(start_date, end_date).select('temperature_2m').mean().subtract(273.15).clip(thailand)
temp_norm = temp.unitScale(20, 30).clamp(0, 1)

# ---------------------------------------------------------
# 3. สร้างฟังก์ชัน: Carrying Capacity (พื้นที่สำหรับ 300 คน)
# ---------------------------------------------------------
def apply_population_rule(suitability_image, min_score=0.6, min_pixels=300):
    habitable = suitability_image.gt(min_score)
    patch_size = habitable.connectedPixelCount(min_pixels, False)
    return suitability_image.updateMask(patch_size.gte(min_pixels))

# ---------------------------------------------------------
# 4. วิเคราะห์เผ่าป่า (Omaticaya: Forest Na'vi)
# ---------------------------------------------------------
landcover = ee.ImageCollection("ESA/WorldCover/v200").first().clip(thailand)
forest_mask = landcover.eq(10)

base_forest = ndvi_norm.multiply(0.4).add(rain_norm.multiply(0.3)).add(temp_norm.multiply(0.2)).add(elev_norm.multiply(0.1))
strict_forest = base_forest.updateMask(forest_mask) 
final_forest = apply_population_rule(strict_forest)

# ---------------------------------------------------------
# 5. วิเคราะห์เผ่าภูเขา (Ash/Fire: Mountain Na'vi)
# ---------------------------------------------------------
slope_norm_easy = slope.unitScale(0, 15).clamp(0, 1)
base_mountain = elev_norm.multiply(0.5).add(slope_norm_easy.multiply(0.5))
strict_mountain = base_mountain.updateMask(elevation.gt(150))
final_mountain = apply_population_rule(strict_mountain, min_score=0.5, min_pixels=200)

# ---------------------------------------------------------
# 6. วิเคราะห์เผ่าน้ำ (Metkayina: Ocean Na'vi)
# ---------------------------------------------------------
ocean_mask = elevation.lte(0)
dist_to_coast = ocean_mask.distance(ee.Kernel.euclidean(40000, 'meters'))
coastal_prox = ee.Image(1).subtract(dist_to_coast.divide(40000)).clamp(0, 1)

base_ocean = coastal_prox.multiply(0.8).add(ee.Image(1).subtract(elev_norm).multiply(0.2))
strict_ocean = base_ocean.updateMask(coastal_prox.gt(0))
final_ocean = apply_population_rule(strict_ocean, min_score=0.5, min_pixels=400)

# ---------------------------------------------------------
# 7. แสดงผลบนแผนที่ (Visualization)
# ---------------------------------------------------------
vis_params = {'min': 0, 'max': 1, 'palette': ['red', 'orange', 'yellow', 'limegreen', 'darkgreen']}

Map.addLayer(final_forest, vis_params, "Omaticaya_Forest", False)
Map.addLayer(final_mountain, vis_params, "AshTribe_Mountain", False)
Map.addLayer(final_ocean, vis_params, "Metkayina_Ocean", True)

empty_style = {'fillColor': '00000000', 'color': 'black', 'width': 1.5}
Map.addLayer(thailand.style(**empty_style), {}, 'Thailand Boundary')

Map
