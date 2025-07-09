from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="vaxlab-report",
    version="0.9",
    author="Chang Lab",
    author_email="",
    description="RNA sequence evaluation and HTML report generation toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ChangLabSNU/VaxLab-report",
    packages=find_packages(),
    package_data={
        "vaxlab_report": [
            "report_template/*.html",
            "report_template/*.css",
            "data/*.py",
            "data/*.json",
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
)