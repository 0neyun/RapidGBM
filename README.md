# RapidGBM

<img align="left" src="https://github.com/user-attachments/assets/05c4d00c-ae6a-44cf-bf41-daa95e63b75b" alt="RapidGBM Logo" width="270"/>

**RapidGBM** is a lightweight, interactive web application and a module of HEtools, designed for rapid preview of Fermi/GBM data. Its core capability is orbit-based visibility calculation: users can instantly determine GBM’s line-of-sight to any sky coordinate at a given time. As soon as continuous TTE data are released, RapidGBM automatically generates detector response files and performs spectral analysis. It is expected to be of assistance during SVOM-BA and EP-TA duty shifts.

Code paper: [Wang et al. 2023](https://arxiv.org/abs/2303.11083).

The tool makes use of publicly available [Fermi GBM tools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/) APIs, and the orbital calculation draw upon the [osv](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/Fermi_GBM_OrbitalBackgroundTool.pdf) from Fermi User Contributions.

<br clear="left"/>

---
