# PRM with an inclusion list in R
## Overview
This R-based tool helps analyzing prm data when collected with an inclusion list. It solves the central problem that DDA-based acquisition using an inclusion list does not "order" the prm MS2 acquisitions in neat transitions like MRM and simply requires selecting MS2 scans of the inclusion list precursor ions from the raw data. The code operates in two phases:
1. **Definition and Selection Phase:** Initially, the tool defines spectra for each target and selects the top 10 precursor ion monitoring (PRM) MS2s based on specific criteria.
2. **Application Phase:** In the subsequent run, it applies these selected PRMs to every analyte and integrates them for comprehensive analysis.

The code is run twice:
- Once to define spectra for each target and select the top 10 PRM MS2s.
- A second time to apply those PRMs to every analyte and integrate them.

## Features
- Automated spectra definition for targeted analysis.
- Selection of top 10 PRM MS2s for focused analysis.
- Integration of selected PRMs across analytes for detailed insights.

## Getting Started

### Prerequisites
- Docker installed on your machine. [Get Docker](https://docs.docker.com/get-docker/)

### Installation and Running with Docker
1. Clone this repository to your local machine.
2. Navigate to the cloned repository directory.
3. Build the Docker image:  ```docker build -t prm . ```
4. Run the tool inside a Docker container: ```docker run --rm prm```

### Without Docker
- If you prefer not to use Docker, dont!

## Contributing
We welcome contributions to enhance the functionality and efficiency of this tool. If you would like to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make your changes and commit them to your branch.
4. Push your branch to your forked repository.
5. Open a pull request from your branch to the main repository.

Please ensure that your code follows our coding standards and includes appropriate tests. We appreciate your contributions!

## Please cite the original paper if used in your work:
Marney, L., Choi, J., Magana, A. A., Yang, L., Techen, N., Alam, N. M., Brandes, M., Soumyanath, A., Stevens, F., & Maier, C. (2024). Liquid Chromatography-Mass Spectrometry Quantification of Phytochemicals in Withania somnifera Using Data-Dependent Acquisition, Multiple-Reaction-Monitoring, and Parallel-Reaction-Monitoring with an Inclusion List. *Front. Chem., Sec. Analytical Chemistry*, 12, 1003. doi: 10.3389/fchem.2024.1373535

## License
This project is licensed under the MIT License - see the `LICENSE` file for more details.

## Contact
For support or to report issues, please file an issue on the GitHub repository or email Luke Marney.
