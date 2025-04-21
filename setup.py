from setuptools import setup, find_packages

setup(
    name='vaxlab-report',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'jinja2',
        'pytrf',
        'tqdm',
        'tabulate',
        'pandas',
        'plotly',
        'pylru'
    ],
    entry_points={
        'console_scripts': [
            'vaxlab-report = vaxlab_report.report_only:main',
            'vaxlab-evaluate = vaxlab_report.evaluate_only:main'
        ]
    },
)
