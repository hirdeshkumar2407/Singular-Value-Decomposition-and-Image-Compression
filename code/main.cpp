#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/IterativeLinearSolvers> 
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <tuple>
#include <utility>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <cmath>
#include <array>
//#include "headers/iterative.h"

#define STB_IMAGE_IMPLEMENTATION
#include "headers/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "headers/stb_image_write.h"
#include "headers/filters.h"

using namespace Eigen;
using namespace std;

// Function to get image dimensions and load the image data
tuple<int, int, int, unsigned char *> getImageDimensions(char *input_image_path)
{
    
    int width, height, channels;
    unsigned char *image_data = stbi_load(input_image_path, &width, &height, &channels, 0);
    if (!image_data)
    {
        cerr << "Error: Could not load image " << input_image_path << endl;
        exit(1);
    }

    return make_tuple(width, height, channels, image_data);
}


// Task1 Matrix and Transpose multiplcation and norm 
MatrixXd task1_productnorm(MatrixXd &A)
{
   MatrixXd matrix_transpose = A.transpose();
    MatrixXd product = matrix_transpose * A;
    int norm = product.norm();

    cout << "\n--------TASK 1----------\n";
    cout << "Norm of A and A transpose: \n" << norm << endl;
    return product;
}

// Task2 Find the two largest singular values
vector<double> task2_findLargestSingularValues(const MatrixXd& A) {
   

    // Step 1: Solve for eigenvalues of A^T * A
    SelfAdjointEigenSolver<MatrixXd> eigenSolver(A);
  

    // Step 2: Get eigenvalues and compute singular values
    VectorXd eigenvalues = eigenSolver.eigenvalues();
    vector<double> singularValues;

    for (int i = 0; i < eigenvalues.size(); ++i) {
        singularValues.push_back(sqrt(eigenvalues[i]));
    }

    // Step 3: Sort singular values in descending order
    sort(singularValues.rbegin(), singularValues.rend());

    // Step 4: Report the two largest singular values
    cout << "\n--------TASK 2----------\n";
    cout << "The two largest singular values by Eigen library method are: "
         << singularValues[0] << " and " << singularValues[1] << endl;
     

     vector<double> result(singularValues.begin(), singularValues.begin() + 2);
    return result;
}

void exportMatrixToMTX(const MatrixXd &A, const std::string &filename)
{
    SparseMatrix<double> ATA_sparse = A.sparseView();
     
    // Export to .mtx file
    if (saveMarket(ATA_sparse, filename))
    {
        std::cout << "A and A transpose matrix successfully exported to " << filename << std::endl;
    }
    else
    {
        std::cerr << "Failed to export A and A transpose matrix." << std::endl;
        ;
    }
}

double getCMDValue(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    
    // Read the command output
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    // Convert the result to double
    std::stringstream ss(result);
    double value;
    ss >> value; // Read the first double from the output string
    
    return value;
}


// Task-3 

// highest eigenvalue of the matrix product by LIS 
double highestEigenValueLIS(MatrixXd A){
 
    exportMatrixToMTX(A, "ATA-product.mtx");
 const char *computelargesteigen ="mpirun -n 1 ./iterativesolver ATA-product.mtx eigvec.txt hist.txt -e pi -tol 1.0e-8 ";
  system(computelargesteigen);

  double largestEigenValue=getCMDValue("mpirun -n 1 ./iterativesolver ATA-product.mtx eigvec.txt hist.txt -e pi -tol 1.0e-8  | grep 'Power: eigenvalue' | awk '{print $4}'");

cout << "The largest eigenvalue of the A and A transpose matrix product by LIS library method is: " << largestEigenValue << endl;

return largestEigenValue;
}

// Compare the result of the highest eigenvalue of the matrix product by LIS and task 2 Single Value

void areValuesClose(double A, double B){
    
      A=sqrt(A);
    cout << "Sqaure Root Eigenvalue by LIS: " << A << endl;
  

    cout << "\nSingle Value by Eigen library: " << B << endl;
    double diff = fabs(A - B);

    if (diff < 0.0001) {
        cout << "The two values are approximately equal." << endl;
    } else {
        cout << "The two values are not equal." <<diff << endl;
    }
}



 

// export image in .png format
void exportimagenotnormalise(unsigned char *image_data, const string &image_name, const string &ext, const MatrixXd &image, int width, int height)
{

    // Create an Eigen matrix to hold the transformed image in unsigned char format
    Matrix<unsigned char, Dynamic, Dynamic, RowMajor> transform_image(height, width);

    // Map the matrix values to the grayscale range [0, 255]
    transform_image = image.unaryExpr([](double val) -> unsigned char
        {
            return static_cast<unsigned char>(std::min(std::max(val, 0.0), 255.0)); // Ensure values stay in [0, 255]
        });

    string output_image_path = "results/" + image_name + "." + ext;

    cout << "image created: " << output_image_path << endl;
    if (stbi_write_png(output_image_path.c_str(), width, height, 1, transform_image.data(), width) == 0)
    {
        cerr << "Error: Could not save grayscale image" << endl;
    }
}

// Function to convert MatrixXd to a VectorXd
VectorXd convertMatrixToVector(const MatrixXd &image_matrix)
{
    VectorXd flattened_image = Map<const VectorXd>(image_matrix.data(), image_matrix.size());
    // Use Map<const VectorXd> instead of Map<VectorXd> because image_matrix.data() returns a const double* (read-only pointer).
    // Eigen allows you to use the Map class to reinterpret the matrix's underlying data as a vector without explicitly copying the data.
    // image_matrix.data() gives a pointer to the matrix's data stored in memory.
    // image_matrix.size() gives the total number of elements in the matrix (i.e., rows * cols).
    // Map<VectorXd> constructs a VectorXd that views the matrix's data as a 1D vector.
    return flattened_image;
}

MatrixXd spMatrixVectorMultiplication(const SparseMatrix<double> &sp_matrix, const VectorXd &vec, int width, int height)
{
    VectorXd result = sp_matrix * vec;
    MatrixXd resultmat = Map<MatrixXd>(result.data(), height, width);
    return resultmat;
}

bool isSymmetric(const SparseMatrix<double>& matrix) {
    if (matrix.rows() != matrix.cols()) {
        return false;  // A non-square matrix cannot be symmetric
    }

    SparseMatrix<double> transposeMatrix = matrix.transpose();
    return matrix.isApprox(transposeMatrix);
}


void exportVectortomtx(const VectorXd &v, const char*  filename)
{
    int n = v.size();
    FILE* out = fopen(filename,"w");
    fprintf(out,"%%%%MatrixMarket vector coordinate real general\n");
    fprintf(out,"%d\n", n);
    for (int i=0; i<n; i++) {
        fprintf(out,"%d %f\n", i ,v(i));
    }
    fclose(out);
        cout << filename << " vector exported successfully"  << endl;

}

VectorXd loadVectorFromMTX(const string &filename) {
    // Load the vector from file
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return VectorXd();
    }

    string line;
    getline(file, line); // Skip the first line (MatrixMarket header)

    int numEntries;
    file >> numEntries;

    VectorXd vec(numEntries);
    for (int i = 0; i < numEntries; ++i) {
        int index;
        double value;

        file >> index >> value;
        vec[index - 1] = value;  // Convert to 0-based indexing
    }

    file.close();
    return vec;
}
 
MatrixXd exportVectorMTXto2DMatrix(const char*  filename, int height, int width)
{
    VectorXd vec=loadVectorFromMTX(filename);
 
    MatrixXd resultmat = Map<MatrixXd>(vec.data(), height, width);
    return resultmat;
}


// Main function
int main(int argc, char *argv[]){

    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <image_path>" << endl;
        return 1;
    }

    char *input_image_path = argv[1];

    // Load image dimensions and data
    auto [width, height, channels, image_data] = getImageDimensions(input_image_path);

    // Prepare Eigen matrix for the image
    MatrixXd image_matrix(height, width);
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            int index = (i * width + j) * channels;
            image_matrix(i, j) = static_cast<double>(image_data[index]); // Keep values in [0, 255]
        }
    }

    // -- Task 1

    MatrixXd product = task1_productnorm(image_matrix);

    
    // -- Task 2
    vector<double> singularValues=task2_findLargestSingularValues(product);
 
    // -- Task 3
    cout << "\n--------TASK 3----------\n";
   double A_Eigen_Value=highestEigenValueLIS(product);

    areValuesClose( A_Eigen_Value,singularValues[0]); 
 
     return 0;
}