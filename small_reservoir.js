// Array containing country information
var countriesInfo = [
  { name: 'Afghanistan', asset: 'afg' },
  { name: 'Algeria', asset: 'alg' },
  { name: 'Bahrain', asset: 'bah' },
  { name: 'Djibouti', asset: 'dji' },
  { name: 'Egypt', asset: 'egy' },
  { name: 'Eritrea', asset: 'eri' },
  { name: 'Ethiopia', asset: 'eth' },
  { name: 'Gaza Strip', asset: 'gzs' },
  { name: 'Iran', asset: 'ira' },
  { name: 'Iraq', asset: 'irq' },
  { name: 'Israel', asset: 'isr' },
  { name: 'Jordan', asset: 'jrd' },
  { name: 'Kuwait', asset: 'kuw' },
  { name: 'Lebanon', asset: 'leb' },
  { name: 'Libya', asset: 'lib' },
  { name: 'Mauritania', asset: 'mau' },
  { name: 'Morocco', asset: 'mor' },
  { name: 'Oman', asset: 'omn' },
  { name: 'Pakistan', asset: 'pak' },
  { name: 'Qatar', asset: 'qtr' },
  { name: 'Saudi Arabia', asset: 'sau' },
  { name: 'Somalia', asset: 'som' },
  { name: 'South Sudan', asset: 'stsd' },
  { name: 'Sudan', asset: 'sud' },
  { name: 'Syria', asset: 'syr' },
  { name: 'Tunisia', asset: 'tun' },
  { name: 'Turkey', asset: 'tur' },
  { name: 'United Arab Emirates',  asset: 'uae' },
  { name: 'Western Sahara', asset: 'wsr' },
  { name: 'West Bank', asset: 'wbk' },
  { name: 'Yemen', asset: 'yem' },
];

// Input reservoir size classes 300 to 1,000, 1,001 to 10,000, 10,001 to 20,000 m^2....., change the short (form) accordingly  
var resSize1 = 300;
var resSize2 = {number:1000, short: '1000'};
var scaleNumber = ee.Number(resSize2.number);
var scale = ee.Number(ee.Algorithms.If(scaleNumber.lte(10000), 10, 100));

// Function to process data for a specific country and year
var processCountryYear = function processCountry(countryInfo, year) {
  
  // Load country borders
  var countryBoundary = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
                          .filter(ee.Filter.eq('country_na', countryInfo.name))
                          .geometry();
  
  //Map.addLayer(countryBoundary, {}, 'Country: ' + countryInfo.name, false);
  
  // Sentinel-2 imagery
  var dataset = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
                  .select('B2', 'B3', 'B4', 'B8')
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .filterDate(year + '-01-01', (parseInt(year, 10) + 1) + '-01-01')
                  .median();
  
  // Load cropland for respective countries from the assets
  var agri = ee.Image("projects/ee-sankeerthg35/assets/" + countryInfo.asset + "Cropland");
  //Map.addLayer(agri.selfMask(), {palette: ['red']}, 'Cropland', false);
  
  // NDWI for water identification
  var ndwi = dataset.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var water = ndwi.updateMask(ndwi.select('NDWI').gt(0.2));
  var waterMask = water.reduce(ee.Reducer.mean()).updateMask(agri);
  var water1 = waterMask.gt(0.2)
                  .selfMask()
                  .rename('water');
  //Map.addLayer(water, {min: 0.0, max: 1.0, palette: ['00FFFF', '0000FF']}, 'Water Mask', false);
  
  // Total yearly Water availability
  //var tWaterAr = water.multiply(ee.Image.pixelArea()).reduceRegion({
  //  reducer: ee.Reducer.sum(),
  //  geometry: countryBoundary,
  //  scale: 30,
  //  crs: 'EPSG:4326',
  //  bestEffort: true,
  //  maxPixels: 1e13,
  //});
  
  // Object-based method for small reservoir identification
  var connectedWater = water1.connectedComponents(ee.Kernel.plus(1), 110);
  var connectedWaterSize = connectedWater.select('labels').connectedPixelCount(110).rename('size');
  var pixelArea = ee.Image.pixelArea();
  var objectArea = connectedWaterSize.multiply(pixelArea);
  var areaMask = objectArea.gte(resSize1).and(objectArea.lte(resSize2.number));
  var resLabel = connectedWater.updateMask(areaMask);
  var smallWaterArea = objectArea.updateMask(areaMask);
  // Map.addLayer(smallWaterArea, {palette: ['yellow']}, 'Small ResArea', false);
  
  // Projection and bounding boxes for 50kmx50km grid boxes are defined
  var proj = ee.Projection("EPSG:4326").atScale(50000);
  // Country bbox
  var grid = countryBoundary.bounds().coveringGrid(proj, 50000);
  
  // Number of small reservoirs within the defined areaMask is mapped over the grid defined
  var densityNumber = resLabel.select("labels").reduceRegions({
    collection: grid,
    reducer: ee.Reducer.countDistinct(),
    scale: scale,
    crs: 'EPSG:4326',
    tileScale: 4
  });
  
  // Export results
  Export.table.toDrive({
    collection: densityNumber,
    description: 'X_' + resSize2.short  + '_' + countryInfo.asset + '_' + year,
    folder: countryInfo.asset + year,
    fileFormat: 'SHP'
  });
  
  // Create Density Map in Region -> Cumulative Area 
  // scale 10 is used for areaMask 300-10k sqm and scale of 100 is used for areaMask 10k-1000k sqm
  var waterPixelArea = areaMask.multiply(ee.Image.pixelArea()); 
  // get sum within each polygon 
  var densityArea = waterPixelArea.reduceRegions({ 
    collection: grid, 
    reducer: ee.Reducer.sum(), 
    scale: scale, 
    crs: 'EPSG:4326', 
    tileScale: 4
  }); 
  
  //export result as shapefile to drive 
  Export.table.toDrive({ 
    collection: densityArea, 
    description: 'X2_' + resSize2.short + '_' + countryInfo.asset + '_' + year, 
    folder: countryInfo.asset + year, 
    fileFormat: 'SHP' 
  });
  
  // Return null to satisfy user-defined method return requirement
  return null;
}

var years = ['2016','2017','2018','2019','2020','2021','2022','2023'];

countriesInfo.forEach(function(countryInfo) {
  years.forEach(function(year) {
    processCountryYear(countryInfo, year);
  });
});
