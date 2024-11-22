# Testing methods for local adaptation

This is a repository for the scripts and data related to the project **"A method for identifying spatially divergent selection considering structured populations."**
This is a project in collaboration with Pierre de Villemereuil, Oscar Gaggiotti and Jerome Goudet.

In this project, we tested three methods:
1. **QstFst**: Based on Whitlock and Guillaume's 2009 approach ([DOI: 10.1534/genetics.108.099812](https://doi.org/10.1534/genetics.108.099812)).
2. **DRIFTSEL**: Karhunen et al 2013 method ([DOI: 10.1111/1755-0998.12111](https://doi.org/10.1111/1755-0998.12111)).
3. **LogAV**: A method we developed in this project.

We tested these methods using simulations of neutrally evolving populations under three population structures:
- 139
- Island Model (IM)
- Stepping Stones (SS)

We also tested our method (LogAV) using populations under selection under two population structures: Stepping Stones and Island Model.

## Repository Structure

### `Scripts/`
- **`Job_Scripts/`**: Contains job scripts used for running simulations and analyses.
- **`Simulation_Scripts/`**: Contains simulation configuration files for quantiNemo.
- **`Method_Testing_Scripts/`**: Contains R scripts for testing each method.
- **`Graph_Scripts/`**: Contains R scripts for creating graphs from results.

### `Data/`
- **`Raw/`**: Contains raw data outputs from quantiNemo simulations.
- **`Processed/`**:
  - **`LAVA/`**: Processed results from the LogAV method.
  - **`Driftsel/`**: Processed results from DRIFTSEL.
  - **`QstFst/`**: Processed results from QstFst.
