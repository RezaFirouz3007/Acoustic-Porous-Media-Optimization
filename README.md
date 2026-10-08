# Acoustic Porous Media Optimization

This repository provides a physics-based computational tool to model and optimize the acoustic absorption of porous materials. It implements the **Johnson-Champoux-Allard (JCA)** and **JCA-Lafarge (JCAL)** macroscopic models, utilizing multi-start L-BFGS-B optimization to inversely identify material parameters from experimental impedance tube data.

## Theoretical Background
The scripts characterize the acoustic behavior of rigid-frame porous media by fitting the following non-acoustic parameters:
* **Porosity** ($\phi$)
* **Flow Resistivity** ($\sigma$)
* **Tortuosity** ($\alpha_\infty$)
* **Viscous Characteristic Length** ($\Lambda$)
* **Thermal Characteristic Length** ($\Lambda'$)
* **Static Thermal Permeability** ($k'_0$) *(JCAL model only)*

## Repository Contents
* `jca_optimization.py`: Fits experimental data using the 5-parameter JCA model.
* `jcal_optimization.py`: Fits experimental data using the 6-parameter JCAL model.
* `jca_jcal_comparison.py`: Comprehensive execution of both models, side-by-side performance comparison, and a full 6-variable sensitivity analysis.
* `my_experiment.xlsx`: Sample experimental dataset containing `frequency` and `absorption` columns. The scripts will automatically write fit metrics, fitted parameters, and plots to new sheets within this workbook.

## Installation & Usage
1. Clone this repository:
   ```bash
   git clone [https://github.com/RezaFirouz3007/Acoustic-Porous-Media-Optimization.git](https://github.com/YOUR-USERNAME/Acoustic-Porous-Media-Optimization.git)