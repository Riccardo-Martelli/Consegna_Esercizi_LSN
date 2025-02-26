# Numerical Simulation Exercises
# Monte Carlo-Based AI and Optimization Algorithms in C++

## Overview
This repository contains a series of **12 projects** focused on **assimilation techniques** using **Monte Carlo methods**, including:
- **Simulated Annealing**
- **Genetic Algorithms**
- **Metropolis Algorithm**

Additionally, I have implemented **basic AI models** using **TensorFlow**, including:
- **Regression Models**
- **Image Recognition** with the **MNIST dataset**

All core implementations are written in **C++**, leveraging **Monte Carlo simulations** to solve various computational tasks efficiently.

## Features
- **Optimization Algorithms**: Implemented Monte Carlo techniques such as **Simulated Annealing** and **Genetic Algorithms**.
- **Metropolis Algorithm**: Used for achieving different computational objectives.
- **Machine Learning with TensorFlow**: Built regressive models and applied deep learning techniques.
- **Image Recognition**: Leveraged the **MNIST dataset** for pattern recognition tasks.
- **Parallel Computing**: Implemented **Parallel Tempering** algorithm using **MPI** for solving the **Traveling Salesman Problem (TSP)**.

## Technologies Used
- **C++**: Core implementation of Monte Carlo-based algorithms.
- **MPI**: Used for parallelizing computations.
- **TensorFlow**: Used for AI models and regression tasks.
- **MNIST Dataset**: Employed for image recognition.
- **Seaborn**: Used for visualization.
- **Glob & ImageIO**: Used for file handling and image processing.

## Installation
To compile and run the C++ projects:
```bash
 g++ -o project_name project_name.cpp -std=c++11
 ./project_name
```

For TensorFlow-based projects:
```bash
pip install tensorflow numpy matplotlib seaborn imageio glob
```

## Usage
Each project is standalone and showcases different aspects of **optimization, AI, and statistical methods**. Explore the subdirectories for detailed documentation on each implementation.

Most of the exercises can be compiled using:
```bash
./execute
```
Otherwise, use:
```bash
make && ./main.cpp
```

Each exercise has a notebook named:
```bash
Exercise_(exercise number).ipynb
```

### Example Exercise: Parallel Tempering Algorithm
One of the exercises requires implementing a **Parallel Tempering algorithm** (Simulated Annealing with multiple temperatures) based on the **Genetic Algorithm** code. The algorithm is parallelized using **MPI** to solve the **Traveling Salesman Problem (TSP)**. Each computing node (up to 10) is assigned a different temperature. The Metropolis algorithm incorporates genetic operators (except crossover) as trial moves, along with a trial move that exchanges paths between adjacent temperatures.

## Future Enhancements
- Extend Monte Carlo methods to real-world applications.
- Improve computational efficiency with parallel processing.

## Contributions
Contributions and suggestions are welcome! Feel free to open an issue or submit a pull request.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
