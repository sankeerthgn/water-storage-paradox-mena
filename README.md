# Water Storage Paradox – MENA (Code & Repro Pack)

Code to (1) identify agricultural reservoirs in MENA with Google Earth Engine,  
(2) compute reservoir evaporation (MATLAB), and  
(3) track atmospheric moisture footprints (Python).

## Repository layout

/gee/
small\_reservoirs.js
/matlab/
evaporation.m
/python/
atmosphere\_tracking.py


## Data (expected paths & names)
- **GEE**: country cropland assets like `projects/ee-<user>/assets/<code>Cropland` (e.g., `afgCropland`).
- **MATLAB**: `Dec_AvgTemp_Regridded.mat`, `MERRA2_regridded_data_<YEAR>.mat`, and `MENA_area_<YEAR>.shp` (e.g., 2016–2023).
- **Python**:  
  - `data_n/Monthly_Avg_EvapLoss.nc` (variables: `lon`, `lat`, `totEvM<month>`),  
  - `data_in/utrack_climatology_0.5_<MM>.nc` (variable: `moisture_flow`).

## Quick start

### 1) Google Earth Engine (GEE) — reservoir identification
- Open the GEE Code Editor → **New** → **Script**.
- Paste `/gee/small_reservoirs.js`.
- Confirm the `countriesInfo` assets (e.g., `projects/ee-<user>/assets/afgCropland`).
- Set `years` and size class `resSize1` / `resSize2`.
- Run to export SHP grids with reservoir **count** and **cumulative area** per 50 × 50 km cell to Google Drive.

### 2) MATLAB — evaporation from reservoirs
- Start MATLAB in `/matlab`.
- Ensure required files are present for the chosen `year` in the script.
- Run:
  ```matlab
  year = 2023;   % edit in the script if needed
  evaporation
``

* Output: `xdailyEvaporation_MENA_4m_<YEAR>.mat` (daily evaporation per grid cell).

### 3) Python — atmospheric moisture tracking

* Python ≥3.9. Install dependencies:

  ```bash
  cd python
  python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
  pip install -r requirements.txt
  ```
* Run for a month (e.g., July):

  ```bash
  python atmosphere_tracking.py 7
  ```
* Output: `output.nc` with variable `forward_footprint` on a 0.5° grid.

## Repro notes

* Coordinate systems: scripts assume regular lat/lon grids (EPSG:4326).
* Large/raster inputs aren’t included.
* For GEE, adjust `scale` if you change reservoir size classes.

````
