// ================================================================================================================================= //
// SURFACE URBAN HEAT ISLAND (SUHI) & UTFVI TEMPORAL ANALYSIS PIPELINE
// Study Area: Special Capital Region of Jakarta (DKI Jakarta), Indonesia
// Observation Period: Dry Season (June - September / JJAS, 2013–2019)
//
// AUTHOR:
// Rafiq Naufal Shidiq (Undergraduate Thesis, Department of Meteorology, Institut Teknologi Bandung - 2020)
// ITB Central Library (Digilib)
//
// METHODOLOGY OVERVIEW:
// 1. Planetary Cloud Screening via Landsat 8 Simple Cloud Score (Sensitivity-tuned threshold <= 22%)
// 2. Fractional Vegetation Cover (FVC) & Emissivity (epsilon) Modeling via NDVI
// 3. Single-Channel Radiative Transfer Formulation for Land Surface Temperature (LST)
// 4. Urban Thermal Field Variance Index (UTFVI) Extraction & Multi-Temporal Statistical Aggregation
//
// NOTE: Script processes multi-temporal composites annually. Date filters can be toggled per study year.
// ================================================================================================================================= //


// ---------------------------------------------------------------------------------------------------------------------------------
// 1. STUDY AREA DEFINITION (DKI JAKARTA)
// ---------------------------------------------------------------------------------------------------------------------------------
// Sub-district (Kelurahan) administrative boundary of DKI Jakarta
var jakarta = ee.FeatureCollection("users/rafiqns/dki_kelurahan");
jakarta = jakarta.geometry();
Map.centerObject(jakarta, 11);
Map.addLayer(jakarta, { color: 'black' }, 'DKI Jakarta Extent');

var aoi = jakarta;

// Standardized inspection bounding box for multi-temporal comparative export
var batas_snip = ee.Geometry.Polygon(
    [[[106.67573593430123, -6.086757742457927],
    [106.67573593430123, -6.376857333770041],
    [106.98403976730904, -6.376857333770041],
    [106.98403976730904, -6.086757742457927]]], null, false);

Map.addLayer(batas_snip, { color: 'black' }, 'Framing Extent');


// ---------------------------------------------------------------------------------------------------------------------------------
// 2. CLOUD SCREENING ALGORITHM (USGS SIMPLE CLOUD SCORE)
// ---------------------------------------------------------------------------------------------------------------------------------
// Cloud-likelihood score threshold set to 22% based on empirical dry-season sensitivity analysis (June 2013)
// Pixels with cloud score > 22% are masked out from further radiative transfer processing
var cloudMask = function (image) {
    var clouds = ee.Algorithms.Landsat.simpleCloudScore(image).select(['cloud']);
    return image.updateMask(clouds.lte(22));
};


// ---------------------------------------------------------------------------------------------------------------------------------
// 3. NORMALIZED DIFFERENCE VEGETATION INDEX (NDVI) FUNCTION
// ---------------------------------------------------------------------------------------------------------------------------------
var addNDVI = function (image) {
    return image.addBands(image.normalizedDifference(['B5', 'B4']).rename('NDVI'));
};


// ---------------------------------------------------------------------------------------------------------------------------------
// 4. INGEST, SPATIAL FILTER, AND MASK LANDSAT 8 TOA COLLECTION
// ---------------------------------------------------------------------------------------------------------------------------------
// Filtered to WRS Path 122, Row 64 (DKI Jakarta footprint) across dry-season months (JJAS)
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterBounds(aoi)
    .filter(ee.Filter.eq('WRS_PATH', 122))
    .filter(ee.Filter.eq('WRS_ROW', 64))
    .filterDate('2013-06-01', '2013-10-01')
    .filter(ee.Filter.date('2013-08-09', '2013-08-10').not())
    //	.filter(ee.Filter.date('2014-07-27', '2014-07-28').not())
    //	.filter(ee.Filter.date('2015-06-12', '2015-06-13').not())
    //	.filter(ee.Filter.date('2015-07-30', '2015-07-31').not())
    //	.filter(ee.Filter.date('2015-08-15', '2015-08-16').not())
    //	.filter(ee.Filter.date('2015-09-16', '2015-09-17').not())
    //	.filter(ee.Filter.date('2016-06-14', '2016-06-15').not())
    //	.filter(ee.Filter.date('2016-07-16', '2016-07-17').not())
    //	.filter(ee.Filter.date('2016-08-01', '2016-08-02').not())
    //	.filter(ee.Filter.date('2016-09-02', '2016-09-03').not())
    //	.filter(ee.Filter.date('2016-09-18', '2016-09-19').not())
    //	.filter(ee.Filter.date('2017-06-01', '2017-06-02').not())
    //	.filter(ee.Filter.date('2017-08-20', '2017-08-21').not()) // Anomalous artifact scene
    //	.filter(ee.Filter.date('2017-09-05', '2017-09-06').not())
    //	.filter(ee.Filter.date('2018-06-20', '2018-06-21').not())
    //	.filter(ee.Filter.date('2018-07-22', '2018-07-23').not())
    //	.filter(ee.Filter.date('2018-08-07', '2018-08-08').not())
    //	.filter(ee.Filter.date('2018-08-23', '2018-08-24').not())
    //	.filter(ee.Filter.date('2018-09-08', '2018-09-09').not())
    //	.filter(ee.Filter.date('2019-06-07', '2019-06-08').not())
    //	.filter(ee.Filter.date('2019-07-09', '2019-07-10').not())
    //	.filter(ee.Filter.date('2019-08-26', '2019-08-27').not())
    //	.filter(ee.Filter.date('2019-09-27', '2019-09-28').not())
    .map(addNDVI)
    .map(cloudMask);

print('Filtered Landsat 8 Collection:', l8);
// Note: Commented dates represent screened scenes with low spatial coverage (>70% cloud cover)


// ---------------------------------------------------------------------------------------------------------------------------------
// 5. SELECT THERMAL (BAND 10) AND NDVI BANDS
// ---------------------------------------------------------------------------------------------------------------------------------
var l8_2 = l8.select(['B10', 'NDVI']);
var l8_2_NDVI = l8.select(['NDVI']);


// ---------------------------------------------------------------------------------------------------------------------------------
// 6. LAND SURFACE TEMPERATURE (LST) RADIATIVE TRANSFER FORMULATION
// ---------------------------------------------------------------------------------------------------------------------------------
// Derives surface emissivity using fractional vegetation cover derived from NDVI,
// then solves the single-channel radiative transfer equation to retrieve LST in degrees Celsius
var toCelsius2_l8 = function (image) {
    var time = image.get('system:time_start');
    var celsius = image.expression(
        '(B10 / (1 + (10.8 * B10 / 14388) * log((0.004 * (((ndvi - 0.2) / 0.3) * ((ndvi - 0.2) / 0.3)) + 0.986)))) - 273.15',
        {
            'ndvi': image.select('NDVI'),
            'B10': image.select('B10')
        })
        .rename("celsius")
        .set('system:time_start', time);
    return celsius;
};


// ---------------------------------------------------------------------------------------------------------------------------------
// 7. MAP LST RETRIEVAL ACROSS COLLECTION
// ---------------------------------------------------------------------------------------------------------------------------------
var l8_2_2 = l8_2.map(toCelsius2_l8);


// ---------------------------------------------------------------------------------------------------------------------------------
// 8. TEMPORAL MEAN LST COMPOSITE
// ---------------------------------------------------------------------------------------------------------------------------------
var l8_2_2_mean_region_image = l8_2_2.reduce(ee.Reducer.mean());


// ---------------------------------------------------------------------------------------------------------------------------------
// 9. RENDER MEAN LST SPATIAL LAYER
// ---------------------------------------------------------------------------------------------------------------------------------
Map.addLayer(l8_2_2_mean_region_image.clip(aoi),
    { palette: ['green', 'yellow', 'orange', 'red'], min: 30, max: 40, scale: 30 },
    'Mean LST (°C)');


// ---------------------------------------------------------------------------------------------------------------------------------
// 10. BUILD INTERACTIVE LST MAP LEGEND (UI PANEL)
// ---------------------------------------------------------------------------------------------------------------------------------
var vizLST = { min: 30, max: 40, palette: ['green', 'yellow', 'orange', 'red'] };

var legendLST = ui.Panel({
    style: {
        position: 'bottom-left',
        padding: '8px 15px'
    }
});

var legendTitleLST = ui.Label({
    value: 'Mean LST (°C)',
    style: {
        fontWeight: 'bold',
        fontSize: '13px',
        margin: '0 0 4px 0',
        padding: '0'
    }
});
legendLST.add(legendTitleLST);

var lonLST = ee.Image.pixelLonLat().select('latitude');
var gradientLST = lonLST.multiply((vizLST.max - vizLST.min) / 100.0).add(vizLST.min);
var legendImageLST = gradientLST.visualize(vizLST);

legendLST.add(ui.Panel({
    widgets: [ui.Label(vizLST['max'] + ' °C')]
}));

legendLST.add(ui.Thumbnail({
    image: legendImageLST,
    params: { bbox: '0,0,10,100', dimensions: '10x200' },
    style: { padding: '1px', position: 'bottom-center' }
}));

legendLST.add(ui.Panel({
    widgets: [ui.Label(vizLST['min'] + ' °C')]
}));

Map.add(legendLST);


// ---------------------------------------------------------------------------------------------------------------------------------
// 11. GENERATE REGIONAL LST TIMESERIES PROFILE CHART
// ---------------------------------------------------------------------------------------------------------------------------------
var COLOR = {
    DKI: '0000ff',
};

var tempTimeSeries = ui.Chart.image.seriesByRegion({
    imageCollection: l8_2_2,
    regions: aoi,
    reducer: ee.Reducer.mean(),
    band: "celsius",
    scale: 30,
    xProperty: 'system:time_start',
    seriesProperty: 'label'
});

tempTimeSeries.setChartType('ScatterChart');
tempTimeSeries.setOptions({
    title: 'DKI Jakarta Land Surface Temperature Timeseries',
    vAxis: { title: 'LST (°C)' },
    hAxis: { title: 'Observation Date' },
    Width: 0.7,
    pointSize: 2,
    series: {
        0: { color: COLOR.DKI }
    }
});

print(tempTimeSeries);


// ---------------------------------------------------------------------------------------------------------------------------------
// 12. URBAN THERMAL FIELD VARIANCE INDEX (UTFVI) PARAMETERIZATION
// ---------------------------------------------------------------------------------------------------------------------------------
var Ts = l8_2_2_mean_region_image;

// Empirical mean urban surface temperature (Tmean) per study year:
var Tmean = 29.07; // 2013 baseline
//var Tmean = 28.10; // 2014
//var Tmean = 30.13; // 2015
//var Tmean = 26.28; // 2016
//var Tmean = 26.94; // 2017
//var Tmean = 31.45; // 2018
//var Tmean = 32.32; // 2019

print('Regional Mean Surface Temperature (Tmean):', Tmean);


// ---------------------------------------------------------------------------------------------------------------------------------
// 13. COMPUTE UTFVI LAYER: UTFVI = (Ts - Tmean) / Ts
// ---------------------------------------------------------------------------------------------------------------------------------
var utfvi = Ts.subtract(Tmean).divide(Ts).rename('UTFVI');

// 14. Render UTFVI Spatial Layer
var utfviParams = { min: -0.05, max: 0.05, palette: ['blue', 'red', 'green'] };
Map.addLayer(utfvi.clip(aoi), utfviParams, 'UTFVI Layer');

// 15. Regional UTFVI Statistical Summary
var utfvi_mean = utfvi.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
});

print('Regional Mean UTFVI:', utfvi_mean);


// ---------------------------------------------------------------------------------------------------------------------------------
// 16. BUILD INTERACTIVE UTFVI MAP LEGEND (UI PANEL)
// ---------------------------------------------------------------------------------------------------------------------------------
var vizUTFVI = { min: -0.05, max: 0.05, palette: ['blue', 'red', 'green'] };

var legendUTFVI = ui.Panel({
    style: {
        position: 'bottom-right',
        padding: '8px 15px'
    }
});

var legendTitleUTFVI = ui.Label({
    value: 'UTFVI Index',
    style: {
        fontWeight: 'bold',
        fontSize: '13px',
        margin: '0 0 4px 0',
        padding: '0'
    }
});
legendUTFVI.add(legendTitleUTFVI);

var lonUTFVI = ee.Image.pixelLonLat().select('latitude');
var gradientUTFVI = lonUTFVI.multiply((vizUTFVI.max - vizUTFVI.min) / 100.0).add(vizUTFVI.min);
var legendImageUTFVI = gradientUTFVI.visualize(vizUTFVI);

legendUTFVI.add(ui.Panel({
    widgets: [ui.Label(vizUTFVI['max'].toString())]
}));

legendUTFVI.add(ui.Thumbnail({
    image: legendImageUTFVI,
    params: { bbox: '0,0,10,100', dimensions: '10x200' },
    style: { padding: '1px', position: 'bottom-center' }
}));

legendUTFVI.add(ui.Panel({
    widgets: [ui.Label(vizUTFVI['min'].toString())]
}));

Map.add(legendUTFVI);


// ---------------------------------------------------------------------------------------------------------------------------------
// 17. OPTIONAL: MUNICIPAL ADMINISTRATIVE BOUNDARIES
// ---------------------------------------------------------------------------------------------------------------------------------
var bag_jakarta = ee.FeatureCollection("users/rafiqns/dki_kota");
Map.addLayer(bag_jakarta, { color: 'white' }, 'Municipal Administrative Boundaries');
