import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="MenDelSIM",
    version="0.0.1",
    author="Tamar V. Av-Shalom",
    author_email="tavshalom@cmmt.ubc.ca",
    description="Mendelian Desease Simulator package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tamario/mendelsim",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)