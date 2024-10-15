# NLA_Challenge2_SVD_IC

# Singular Value Decomposition and Image Compression

## Project Overview

This project demonstrates the implementation of image compression and noise reduction using **Singular Value Decomposition (SVD)** with **C++** and the **Eigen** library. The focus is on applying SVD to perform tasks such as **image compression**, **noise reduction**, and **eigenvalue computation**, while also exploring the use of the **LIS library** for iterative solvers.

---

## 1. Singular Value Decomposition (SVD)

Singular Value Decomposition is a fundamental technique in linear algebra that factorizes a matrix into three matrices: U, Σ, and V^T. For any m × n matrix A:

A = UΣV^T

Where:
- U is an m × m orthogonal matrix
- Σ is an m × n rectangular diagonal matrix with nonnegative real numbers on the diagonal
- V^T is the transpose of an n × n orthogonal matrix

SVD has numerous applications, including image compression and noise reduction.

---

## 2. Image Compression Using SVD

SVD can be used for image compression by approximating the original matrix with a lower-rank matrix. This is done by keeping only the k largest singular values and their corresponding singular vectors:

A ≈ σ₁u₁v₁^T + σ₂u₂v₂^T + ... + σₖuₖvₖ^T

Where σᵢ are singular values, uᵢ are left singular vectors, and vᵢ are right singular vectors.

### Example: 4x4 Image Matrix

Let's walk through a simple example of SVD-based image compression using a 4x4 grayscale image matrix.

#### Step 1: Original Image Matrix A
```
A = [100  120  130  115]
    [105  125  135  120]
    [110  130  140  125]
    [115  135  145  130]
```

#### Step 2: Compute SVD
After computing SVD, we get U, Σ, and V^T matrices.

#### Step 3: Truncate SVD
For compression, we keep only the first k singular values and vectors. Let's say k = 2.

#### Step 4: Reconstruct Compressed Image
The compressed image is reconstructed using only the first two components:

A' ≈ σ₁u₁v₁^T + σ₂u₂v₂^T

#### Step 5: Resulting Compressed Image Matrix
```
A' ≈ [99   119  131  116]
     [104  124  136  121]
     [110  130  141  126]
     [115  135  146  131]
```

This compressed version approximates the original image while using less data.

---

## 3. Why Numerical Linear Algebra?

SVD and image processing involve complex matrix operations, making efficient numerical methods crucial. **Numerical linear algebra** provides the tools to perform these operations effectively. In this project, we use:

- **Eigen Library**: For matrix operations and SVD computation
- **LIS Library**: For iterative solvers in eigenvalue problems

---

## 4. C++ Implementation with Eigen and LIS

We use the **Eigen** and **LIS** libraries for matrix manipulation, SVD computation, and solving linear systems. Here's an overview of the key libraries and headers used:

```cpp
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <fstream>

// Image processing libraries
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// LIS library headers
#include "lis.h"
```

---

## 5. Task Descriptions

This project implements several tasks related to SVD and image processing:

### Task 1: Loading and Processing the Image
- Load the 'einstein.jpg' image as an Eigen matrix A
- Compute A^T * A and report its Euclidean norm

### Task 2: Eigenvalue Computation with Eigen
- Solve the eigenvalue problem A^T * A * x = λx using Eigen
- Report the two largest computed singular values

### Task 3: Eigenvalue Computation with LIS
- Export A^T * A in matrix market format
- Use LIS to compute the largest eigenvalue with a tolerance of 10^-8
- Compare results with Eigen's computation

### Task 4: SVD Computation
- Perform SVD on matrix A using Eigen
- Report the Euclidean norm of the diagonal matrix Σ

### Task 5: Image Compression
- Compute matrices C and D for k = 40 and k = 80
- Generate compressed images and export as PNG

### Task 6: Checkerboard Image Creation and Noise Addition
- Create a 200x200 black and white checkerboard image
- Add random noise in the range [-50, 50] to each pixel

### Task 7: Noise Reduction using SVD
- Perform SVD on the noisy checkerboard image
- Create compressed versions with k = 5 and k = 10
- Compare original, noisy, and denoised images

---

## 6. Results and Analysis

(After completing the project, summarize key findings here, including:)
- Effectiveness of SVD in image compression
- Comparison of eigenvalue computation methods (Eigen vs. LIS)
- Analysis of noise reduction results
- Visual comparisons of original, compressed, and denoised images

---

## Conclusion

This project demonstrates the power of Singular Value Decomposition in image processing, particularly for compression and noise reduction. It showcases the practical application of numerical linear algebra techniques and the use of efficient libraries like Eigen and LIS for complex matrix operations.

--
## Acknowledgements

This was challenge 2 assigned as part of the course **Numerical Linear Algebra in High-Performance Computing (2024/25)** at **Politecnico di Milano**. We extend our sincere gratitude to:

- **Professor [P. F. Antonietti](https://www.linkedin.com/in/paolaantonietti/?lipi=urn%3Ali%3Apage%3Ad_flagship3_search_srp_all%3BtoYfzDyNQUuaYhVlXkVXMQ%3D%3D)**, for providing excellent guidance throughout the course.
- **Teaching Assistant [Dott. M. Botti](https://www.linkedin.com/in/michele-botti-4707a62a2/?lipi=urn%3Ali%3Apage%3Ad_flagship3_search_srp_all%3BFvI80B0lRXiNyhRyRoR13Q%3D%3D)**, for their support and valuable feedback on this project.

This project has significantly enhanced our understanding of **numerical methods**, **image processing**, and their applications in **high-performance computing**.




## Contributors 
[Hirdesh Kumar](https://www.linkedin.com/in/hirdeshkumar2407/)

[Nadah Khaled](https://www.linkedin.com/in/nadahkhaledd10/)

[Milica Sanjevic](https://www.linkedin.com/in/milica-sanjevic-321392327/)





