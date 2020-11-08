
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CSHAL",
    version="1.0",
    author="Pradeep Varathan Pugalenthi",
    author_email="pradeepvarathanpugalenthi@gmail.com",
    description="A pipeline package to process summary statistics file for DeepSEA analysis",
    long_description="Using the effect and alternate allele, along with the C-SHAL's smooth takeover to convert them to sequences and feed into the DeepSEA model.",
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)