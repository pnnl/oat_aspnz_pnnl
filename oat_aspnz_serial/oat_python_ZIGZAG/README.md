# Open Applied Topology in Python

Open Applied Topology in Python (OAT-Python) is a Python package for fast, user-friendly algebra and topology.  It provides wrappers for the Rust library Open Applied Topology.

## Documentation

The best ways to learn about functions and objects are currently (i) the example notebooks, and (ii) Python's `help()` function.

## Installation

1. Download and install the most recent version of [Rust](https://www.rust-lang.org/).  Make sure your installation is up to date by running `rustup update` in a command shell.

2. Obtain a copy of OAT suite, which contains folders for OAT-Rust and OAT-Python.  Don't move the contents out of this folder.

3. Create a virtual Python environment, e.g. using Anaconda, or open one that you already have.  In this example, let's assume the environment name is `myenv`.  Activate `myenv`, and run

    ```bash
    pip install maturin
    ```

    A number of warning messages may appear; this is normal. 

    **If you already have maturin installed, make sure it is up to date!**

4. With `myenv` activated, CD into the `oat_python` folder, and run

    ```bash
    maturin develop --release
    ```
    
5. OAT-Python should now be installed.  Try running the Jupyter notebooks with `myenv`!
