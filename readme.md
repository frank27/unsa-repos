# Matrix Multiplication: Standard vs. DGEMM routine

This project compares the performance of the **Standard ()** matrix multiplication method against **DGEMM routine **. It features a high-performance C++ backend managed by CMake and a Python-based data visualization suite.

## 📋 Prerequisites

Ensure you have the following installed on your system:

### System Tools

* **C++ Compiler:** GCC, Clang, or MSVC (supporting C++11 or higher).
* **CMake:** Version 3.10 or higher.
* **Make:** For Unix-based environments.

### Python Environment

The visualization script requires **Python 3** and the following libraries:

* `pandas`
* `matplotlib`

---

## 🚀 How to Build and Run

I have provided a Bash script to automate the compilation and execution process. This script handles CMake configuration, building, and automatically locates the generated executable.

### Using the Automation Script

1. **Give the script execution permissions:**
```bash
chmod +x build.sh

```


2. **Run the script** by passing the project directory (usually the current directory) as an argument:
```bash
./build.sh .

```



---

## 📊 Data Visualization

The project includes a Python script, `graph.py`, to visualize the execution time and efficiency data output from the C++ application.

### Usage

Once the C++ executable has finished running and generated a `.csv` data file, generate your charts by running:

```bash
python graph.py

```

> **Note:** Ensure the generated CSV file is in the same directory as `graph.py` or update the script path accordingly.