# Jakarta Surface Urban Heat Island (SUHI) Analysis using Google Earth Engine

[![GEE](https://img.shields.io/badge/Platform-Google%20Earth%20Engine-green.svg)](https://earthengine.google.com/)
[![Satellite](https://img.shields.io/badge/Satellite-Landsat%208%20OLI%2FTIRS-blue.svg)](https://landsat.gsfc.nasa.gov/satellites/landsat-8/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Official Google Earth Engine JavaScript pipeline accompanying the undergraduate honours thesis:
> **"Temporal Analysis of the Surface Urban Heat Island Phenomenon in DKI Jakarta Using Google Earth Engine"**  
> *(Analisis Fenomena Surface Urban Heat Island secara Temporal di DKI Jakarta Menggunakan Google Earth Engine)*  
> **Author:** Rafiq Naufal Shidiq (Department of Meteorology, Institut Teknologi Bandung / ITB, 2020)  
> **Thesis Record:** [ITB Central Library Digilib](https://digilib.itb.ac.id/gdl/view_data/analisis-fenomena-surface-urban-heat-island-secara-temporal-di-dki-jakarta-menggunakan-google-earth-engine/rafiq-naufal-shidiq)  
> **Interactive GEE Code Editor Script:** [Open in Google Earth Engine](https://code.earthengine.google.com/27c8a7e3173844e1f8a1eb334e9cfcb4)

---

## 📌 Research Overview

This study investigates the multi-year spatiotemporal dynamics of the **Surface Urban Heat Island (SUHI)** phenomenon across the special capital region of Jakarta (**DKI Jakarta**, covering Central, North, West, South, and East Jakarta) during the dry season (June–September / JJAS) across a 7-year observational window (2013–2019).

### Methodology Summary:
1. **Satellite Data Acquisition:**
   - Evaluated 52 Landsat 8 OLI/TIRS satellite scenes (WRS-2 Path 122, Row 64).
   - Applied **Simple Cloud Score Algorithm** (threshold $\le 22\%$), selecting 30 representative scenes with $> 30\%$ valid terrestrial pixels.
2. **Land Surface Temperature (LST) Derivation:**
   - Converted Top of Atmosphere (TOA) spectral radiance from TIRS Band 10 ($10.60 - 11.19\,\mu\text{m}$) into at-satellite Brightness Temperature ($T_b$) via Planck's radiative transfer inversion:
     $$T_b = \frac{K_2}{\ln\left(\frac{K_1}{L_\lambda} + 1\right)}$$
   - Derived surface emissivity ($\varepsilon$) dynamically utilizing the Sobrino vegetation proportion model ($P_v$ from NDVI):
     $$\varepsilon = m \cdot P_v + n$$
   - Retrieved LST in degrees Celsius ($T_s$) using the Artis & Carnahan atmospheric correction:
     $$T_s = \frac{T_b}{1 + \left(\frac{\lambda \cdot T_b}{\rho}\right) \ln(\varepsilon)} - 273.15$$
3. **Urban Thermal Field Variance Index (UTFVI):**
   - Quantified ecological thermal comfort and environmental stress zones across five administrative municipalities:
     $$\text{UTFVI} = \frac{T_s - T_{\text{mean}}}{T_s}$$
   - Stratified into 6 ecological comfort classes: None ($<0$), Weak ($0 - 0.005$), Moderate ($0.005 - 0.010$), Strong ($0.010 - 0.015$), Stronger ($0.015 - 0.020$), and Strongest ($>0.020$).
4. **Climate Teleconnection Correlation:**
   - Correlated annual dry-season LST anomalies with the **ENSO Oceanic Niño Index (Niño 3.4)**, demonstrating a positive teleconnection coefficient of $r = 0.41$.

---

## 🚀 How to Run in Google Earth Engine

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Create a new script and copy-paste the contents of [`jakarta_suhi_gee.js`](jakarta_suhi_gee.js).
3. Click **Run** to execute the pipeline. The interactive map will render:
   - Annual multi-temporal dry-season LST composites (2013–2019).
   - Dynamic 6-class UTFVI ecological comfort stratification.
   - Interactive inspector and annual zonal thermal charts.

---

## 📁 Repository Structure

```
jakarta-suhi-gee/
├── jakarta_suhi_gee.js    # Complete GEE JavaScript pipeline
├── README.md              # Research documentation & radiative transfer equations
├── .gitignore             # Git ignore rules
└── LICENSE                # MIT License
```

---

## 📄 Citation

If you utilize this methodology or script in your academic work, please cite:

```bibtex
@bachelorsthesis{shidiq2020suhi,
  author       = {Shidiq, Rafiq Naufal},
  title        = {Temporal Analysis of the Surface Urban Heat Island Phenomenon in DKI Jakarta Using Google Earth Engine},
  school       = {Department of Meteorology, Faculty of Earth Sciences and Technology, Institut Teknologi Bandung},
  year         = {2020},
  address      = {Bandung, Indonesia},
  note         = {ITB Central Library Digilib}
}
```

---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
