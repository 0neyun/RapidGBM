# RapidGBM

<img align="left" src="https://github.com/user-attachments/assets/ed371fb2-db1f-4c36-b4c3-e2fe691ec306" alt="RapidGBM Logo" width="270"/>

RapidGBM is a lightweight, interactive web application and a module of HEtools, designed for rapid preview of Fermi/GBM data. Its core capability is orbit-based visibility calculation: users can instantly determine GBM’s line-of-sight to any sky coordinate at a given time. As soon as continuous TTE data are released, RapidGBM automatically generates detector response files and performs spectral analysis. It is expected to be of assistance during SVOM-BA and EP-TA duty shifts.

Code paper: [Wang et al. 2023](https://arxiv.org/abs/2303.11083).

The tool leverages existing APIs from the [Fermi GBM tools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/), and the orbital calculation is based on the [osv](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/Fermi_GBM_OrbitalBackgroundTool.pdf) from Fermi User Contributions.

<br clear="left"/>

---
