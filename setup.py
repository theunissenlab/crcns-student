from setuptools import setup, find_packages


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name = "crcns-student",
    version = "0.1.2",
    packages = find_packages(exclude=["examples*", "test*"]),
    scripts = ["manage.py"],
    include_package_data = True,
    zip_safe = False,
    description = "computational neuroscience workshop series tutorial tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author = "Frederic Theunissen",
    url = "https://github.com/theunissenlab/crcns-student",
    classifiers = [
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.0",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        "download",
        "nbgrader",
        "matplotlib",
        "pandas",
        "PyQt5",
        "tqdm",
        "scipy",
    ],
    entry_points="""
        [console_scripts]
        sep=manage:cli
    """,
)
