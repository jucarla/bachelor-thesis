## Decoupling of Temperature and Strain in Fiber Optics-Based Sensors Integrated into Battery Cells

### Repository Overview

This repository contains the most relevant Python scripts and supporting files developed for my Bachelor Thesis, titled **"Decoupling of Temperature and Strain in Fiber Optics-Based Sensors Integrated into Battery Cells."** The thesis was presented as a requirement for obtaining the degree of **Bachelor in Control and Automation Engineering** from **Centro Universitário SENAI CIMATEC - Salvador, Bahia, Brazil**. All the experiments were conducted in collaboration with the Fraunhofer Institute for Production Technology IPT and focuses on developing advanced sensing methodologies for lithium-ion battery cells.

### Project Summary

Lithium-ion batteries play a critical role in enabling sustainable energy solutions. Monitoring their internal conditions, such as temperature and strain, is essential for improving safety, longevity, and performance. This project explores the application of fiber optics-based sensors and a decoupling algorithm to measure temperature and strain independently.

### Repository Structure

```
bachelor-thesis/
│
├── scripts/
│   ├── data_acquisition.py     # Scripts for collecting data from sensors
│   ├── spectral_analysis.py    # Tools for processing spectral shift data
│   ├── decoupling_algorithm.py # Implementation of the decoupling algorithm
│   ├── visualization.py        # Tools for data visualization
│   └── utilities.py            # Common helper functions
│
├── data/
│   ├── raw/                    # Raw data collected from experiments
│   ├── processed/              # Processed data used for analysis
│   └── results/                # Output of analysis and simulations
│
├── notebooks/
│   ├── kt_determination.ipynb  # Jupyter notebook for KT coefficient analysis
│   └── battery_test.ipynb      # Notebook for validating decoupling in battery cells
│
├── figures/
│   ├── spectral_shift.png      # Example plots of spectral shift vs temperature
│   └── setup_diagram.png       # Diagram of experimental setup
│
└── README.md                   # This README file
```

### Key Features

- **Temperature Sensitivity Analysis:** Scripts for determining the temperature sensitivity coefficient (KT) of boron and germanium-doped fibers.
- **Decoupling Algorithm:** An algorithm to separate strain and temperature effects using spectral data.
- **Data Processing Pipelines:** Automated tools for filtering and analyzing raw data collected during experiments.
- **Visualization Tools:** Scripts to generate plots and visuals for thesis documentation and presentations.

### Requirements

- Python 3.8 or higher
- Required libraries:
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `pandas`
  - `seaborn`
  - `jupyter`


### Thesis Citation

If you use any of the scripts or data from this repository, please cite my thesis:

**Juliana Carla Santos da Silva.** *Decoupling of Temperature and Strain in Fiber Optics-Based Sensors Integrated into Battery Cells.* Centro Universitário SENAI CIMATEC, November 2024.

### Acknowledgments

This project was conducted in collaboration with:
- Fraunhofer Institute for Production Technology IPT
- Supervisors: Caroline Girmen and Paulo Andrade Souza 

### License

This repository is licensed under the [MIT License](LICENSE).
