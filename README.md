# Black-Scholes Option Pricing Project

This repository contains a C++ implementation of a **Black-Scholes Partial Differential Equation (BS PDE) model** for option pricing, along with a **VBA Excel macro-enabled workbook** for enhanced interaction and visualization.

## Overview

The project aims to solve the **BS PDE** numerically to price both **European and American options** and compute the **Greeks** (Delta, Gamma, Theta, Vega, and Rho). The solution integrates a **C++ backend** for the core numerical computations with a **VBA-based Excel interface**, making it accessible and interactive.

### Black-Scholes Partial Differential Equation (BS PDE)

The Black-Scholes PDE is a fundamental equation for option pricing, given by:

<img width="275" alt="image" src="https://github.com/user-attachments/assets/15d87c0b-3d7c-4ec8-9ba6-1591ff26497d" />


Where:
- \( V(S, t) \) is the option price as a function of the underlying asset price \( S \) and time \( t \).
- \( sigma \) is the volatility of the underlying asset.
- \( r \) is the risk-free interest rate.

For a **European call option**, the boundary condition at maturity \( t = T \) is:

\[
V(S, T) = max(S - K, 0)
\]

And for a **European put option**:

\[
V(S, T) = max(K - S, 0)
\]

This PDE is solved numerically in the C++ code using a **Crank-Nicholson finite difference method**, discretizing both time and space to obtain the option prices and Greeks.

### Mathematical Approach

The solution of the **BS PDE** involves several key numerical methods and mathematical steps:

1. **Discretization of the PDE**:
   - We divide the time interval \([0, T]\) into \(M\) time steps, and the asset price range \([0, S_{\text{max}}]\) into \(N\) discrete points.
   - Using the **Crank-Nicholson method**, we approximate the derivatives in the PDE with finite differences, leading to a system of linear equations at each time step.

2. **Tridiagonal Matrix Algorithm (Thomas Algorithm)**:
   - The discretization leads to a tridiagonal matrix system, which we solve efficiently using the **Thomas algorithm**.
   - This algorithm reduces the computational complexity to \(O(N)\) per time step, making the numerical solution practical.

3. **Boundary Conditions**:
   - For European options, we apply the terminal condition at maturity \(t = T\) and solve backward in time.
   - For American options, we account for the early exercise feature by adjusting the solution to satisfy the early exercise condition at each step.

4. **Calculation of Greeks**:
   - The **Greeks** are calculated using finite difference approximations:
   <img width="407" alt="image" src="https://github.com/user-attachments/assets/06538df8-6890-4bd1-a402-2a09e0513232" />

   - Additional Greeks, such as **Vega** and **Rho**, are obtained by perturbing the volatility and interest rate, respectively, and recalculating the option price.

This mathematical approach ensures that the numerical solution is stable and convergent, providing accurate option prices and Greeks.




### C++ Implementation (`code.cpp`)

The main C++ file implements a `BlackScholesSolver` class, which:
- **Discretizes the BS PDE using the Crank-Nicholson method** for time-space grids.
- Solves the **tridiagonal matrix system using the Thomas algorithm**.
- **Calculates option prices** for European and American options.
- **Computes the Greeks** to assess the sensitivity of the option prices.
- **Provides a comparison with the analytical Black-Scholes solution** for validation.

The solver handles different types of options:
- **European Options (Call/Put)**
- **American Options (Call/Put)**

The program can be executed from the command line with the following parameters:
```
Usage: test3.exe N M S0 T t K sigma r IsCall IsAmerican
```
Where:
- `N`: Number of space points
- `M`: Number of time points
- `S0`: Initial asset price
- `T`: Time to maturity
- `t`: Current time
- `K`: Strike price
- `sigma`: Volatility
- `r`: Risk-free rate
- `IsCall`: 1 for Call option, 0 for Put option
- `IsAmerican`: 1 for American option, 0 for European option

### VBA Integration (`pricing_option.xlsm`)

The Excel file includes VBA macros that:
- **Take user input parameters** for option pricing.
- **Execute the C++ program through the command line**.
- **Import the results back into Excel**, displaying the option prices and Greeks over a range of underlying asset prices.

This interaction between **Excel (VBA)** and **C++** creates a practical and efficient solution for users to run complex **BS PDE** computations with a familiar and intuitive interface.

### How It Works

1. **Users input the option parameters in Excel.**
2. The **VBA macro launches the C++ executable** with the provided parameters.
3. The **C++ program solves the BS PDE** and outputs the results.
4. **Results are automatically imported into Excel** and displayed in a user-friendly format.

## Results
<img width="575" alt="image" src="https://github.com/user-attachments/assets/8b3e7e6e-81a4-4250-99c5-4cce26f41c06" />


In the input we enter all the parameters of the option we are interested in.



Then the C++ code is executed and provides the price of the option with it's greeks : 
<img width="576" alt="image" src="https://github.com/user-attachments/assets/54aa2fea-c1b9-4c71-9f06-fdcdec2748b7" />



In the output we also have some charts displaying the price of the option with respect to the value of the underlying and also the value of the delta with respect to the value of the underlying, both of them are created thanks to a table of points imported into Excel after the C++ code is executed 



<img width="359" alt="image" src="https://github.com/user-attachments/assets/a71067fc-6529-4943-8404-512f7f97b37e" />

<img width="351" alt="image" src="https://github.com/user-attachments/assets/6a0881f8-b5c6-42ee-8b1f-7ff2e436d68b" />





