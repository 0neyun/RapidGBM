# RapidGBM

<img align="left" src="https://github.com/user-attachments/assets/05c4d00c-ae6a-44cf-bf41-daa95e63b75b" alt="RapidGBM Logo" width="270"/>

**RapidGBM** is a lightweight, interactive web application and a module of HEtools, designed for rapid preview of Fermi/GBM data. Its core capability is orbit-based visibility calculation: users can instantly determine GBM’s line-of-sight to any sky coordinate at a given time. As soon as continuous TTE data are released, RapidGBM automatically generates detector response files and performs spectral analysis. It is expected to be of assistance during SVOM-BA and EP-TA duty shifts.

Code paper: [Wang et al. 2025](https://arxiv.org/abs/2506.20532).

The tool makes use of publicly available [Fermi GBM tools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/) APIs, and the orbital calculation draw upon the [osv](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/Fermi_GBM_OrbitalBackgroundTool.pdf) from Fermi User Contributions. Response functions are generated using the [GBM Response Generator](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/DOCUMENTATION.html) provided by the Fermi team. Spectral fitting is performed using [PyXspec](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html). The interactive web interface is built with the easily deployable [Streamlit](https://streamlit.io) framework.


We provide a publicly accessible web application that includes this function: **[hetools.xyz](http:hetools.xyz)**



<br clear="left"/>

---
