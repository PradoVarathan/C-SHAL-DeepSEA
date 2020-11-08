

![LinkedIn][linkedin-shield]


<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h3 align="center">C-SHAL : Deep SEA Pipeline Tool</h3>

  <p align="center">
    Chromatin-Sequence Hot encoded Analysis PileLine with DeepSEA network
    <br />
    <a href="https://github.com/PradoVarathan/C-SHAL-DeepSEA"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/PradoVarathan/C-SHAL-DeepSEA">View Demo</a>
    ·
    <a href="https://github.com/PradoVarathan/C-SHAL-DeepSEA/issues">Report Bug</a>
    ·
    <a href="https://github.com/PradoVarathan/C-SHAL-DeepSEA/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
* [Usage](#usage)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project
  This pipeline is built for the purpose of analysis directly from a summary statistics file for a Genome Wide Analysis Study with the tab seperated format containing columns with labels having :
* SNP - single nucleotide variants as either rsids or similar labels
* CHR - chromosome of the SNP
* BP - base pair position of the SNP
* A1 - effect allele of the SNP
* A2 - alternate allele of the SNP
  Using the DeepSEA model, built with the convolution networks in deepleaninig using pytorch and has to be worked with GPU and download the pickle file for the best hyperparameters that can be found at [GoogleLink](https://drive.google.com/file/d/1ogbiyyAc32BNDDTAao8INIJeVWnnube-/view?usp=sharing). Make sure to have this file in the same directory you run the whole script.

### Built With

* DeepLearning Packages - Pytorch, Sklearn
* Data Handling packages - Pandas, Numpy, os, Pickle
* CLI Packages - Click
* BioPython



<!-- GETTING STARTED -->
## Getting Started

Clone the whole git repository into your own system, unix-based OS is required for now. Install the required packages using pip and download the best hyperparameters model in the same folder where you run the script.

<!-- USAGE EXAMPLES -->
## Usage

Just using the simple python command with the other mentioned parameters, as follows,
(`python CSHAL.py --ss='Path_To_Summary_Statistics_File' --w=500 --email='ENtrez_account_email_id' --ak='API KEY' --det=choose among log/all/both`)
with that, you should be prompted with multiple options if you did not fill it up right.
* ss --> A summary statistics tab seperated file with columns and headers as SNP,CHR,BP,A1,A2. The columns containing Single Nulceotide Polymorphisms, chormosome, base pair locations, primary allele and secondary allele.
* --w --> The upstream and downstream width from base pair.
* --email --> Email for the Entrez ID to obtain sequences.
* --ak --> API Key from NCBI for faster processing.
* --det --> Options for the detail in the output file. log only gives basic log terms;all provides all 919 labels and values; both provides both the files.

For more detailed explanation on command line, use (`python CSHAL.py --help`)

<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/PradoVarathan/C-SHAL-DeepSEA/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

 PradeepVarathan - [@VarathanPradeep](https://twitter.com/VarathanPradeep) - pradeepvarathanpugalenthi@gmail.com

Project Link: [https://github.com/PradoVarathan/C-SHAL-DeepSEA](https://github.com/PradoVarathan/C-SHAL-DeepSEA)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* [DeepSEA - Zhou, J., Troyanskaya, O. Predicting effects of noncoding variants with deep learning–based sequence model. Nat Methods 12, 931–934 (2015)](https://doi.org/10.1038/nmeth.3547)
* [ReadME Template](https://github.com/othneildrew/Best-README-Template)







<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=flat-square&logo=linkedin&colorB=555
[linkedin-url]: [https://www.linkedin.com/in/pradeep-varathan-pugalenthi-898b31140]
