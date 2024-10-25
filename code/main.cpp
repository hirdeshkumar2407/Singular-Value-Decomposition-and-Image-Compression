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

// export image in .png format
void normalizeExportImage(unsigned char *image_data, const string &image_name, const string &ext, const MatrixXd &image, int width, int height)
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

  double largestEigenValue=getCMDValue("mpirun -n 1  ./iterativesolver ATA-product.mtx eigvec.txt hist.txt -e pi -tol 1.0e-8  | grep 'Power: eigenvalue' | awk '{print $4}'");

cout << "The largest eigenvalue of the A and A transpose matrix product by LIS library method is: " << largestEigenValue << endl;

return largestEigenValue;
}

// Compare the result of the highest eigenvalue of the matrix product by LIS and task 2 Single Value

void areValuesClose(double A, double B){
   double midpoint = (A / 2) ;
      A=sqrt(A);
    cout << "Sqaure Root Eigenvalue by LIS: " << A << endl;
   

    cout << "\nSingle Value by Eigen library: " << B << endl;
    double diff = fabs(A - B);

    if (diff < 0.0001) {
        cout << "The two values are approximately equal." << "Mid-point:" << midpoint << endl;
    } else {
        cout << "The two values are not equal." <<diff << endl;
    }
}

// Task 4 

void task4_shiftIterativeSolver(){
     const char *computelargesteigen ="mpirun -n 1  ./iterativesolver ATA-product.mtx eigvec.txt hist.txt -e ii -tol 1.0e-8 -shift 1.045818e+09" ;
  system(computelargesteigen);
}

JacobiSVD<MatrixXd> SVD(MatrixXd A){
    JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    // Step 1: Get the singular values
    VectorXd singularValues = svd.singularValues();
    
    // Step 2: Compute the Euclidean norm of the singular values
    double euclideanNorm = singularValues.norm();
    
    // Step 3: Print the result
    cout << "The Euclidean norm of the singular values is: " << euclideanNorm << std::endl;

    return svd;

}

// Task 6
pair<MatrixXd, MatrixXd> createCD(JacobiSVD<MatrixXd> svd, int k){
    cout << " -- starting construction of C and D matrices with k= " << k << "--\n";
        MatrixXd U = svd.matrixU();
        MatrixXd V = svd.matrixV();
        VectorXd S = svd.singularValues();

        int rank = S.size();  // This is the rank of the matrix
       // Adjust k if it exceeds the rank
        if (k > rank) {
            cout << "k is greater than the rank of the matrix. Adjusting k to rank: " << rank << endl;
            k = rank;
        }

        cout << "S dimensions: " << S.rows() << " x " << S.cols() << endl;
        MatrixXd c = U.leftCols(k);
        MatrixXd d = V.leftCols(k) * S.head(k).asDiagonal();

        int nonZeroElementsC = c.nonZeros();
        int nonZeroElementsD = d.nonZeros();

        cout << "C dimensions: " << c.rows() << " x " << c.cols() << endl;
        cout << "D dimensions: " << d.rows() << " x " << d.cols() << endl;

        // Print the result
        cout << "Non-zero elements of C with k= " << k << "  are: \n" << nonZeroElementsC << endl;
        cout << "Non-zero elements of D with k= " << k << "  are: \n" << nonZeroElementsD << endl;
   
   return make_pair(c, d);
}

MatrixXd computeExportCompressedImage(pair<MatrixXd, MatrixXd> cd, const string &filename)
{
    MatrixXd C = cd.first;
    MatrixXd D = cd.second;
   
    MatrixXd DT= D.transpose();

    
    cout << "C dimensions: " << C.rows() << " x " << C.cols() << endl;
    cout << "D dimensions: " << DT.rows() << " x " << DT.cols() << endl;
    MatrixXd compressedA = C * DT;

    // cout << "Compressed A: \n" <<compressedA << endl;
    cout << "Compressed Image dimensions: " << compressedA.rows() << " x " << compressedA.cols() << endl;

    // Export the compressed matrix to .mtx file
    normalizeExportImage(nullptr, filename, "png", compressedA, compressedA.cols(),compressedA.rows());

    return compressedA;
}
MatrixXd createCheckerBoard(int img_size,int square_size){

   
 
 MatrixXd checkerboard(img_size, img_size);
    for (int i = 0; i < img_size; ++i)
    {
        for (int j = 0; j < img_size; ++j)
        {
            if (((i / square_size) % 2 == 0 && (j / square_size) % 2 == 0) || ((i / square_size) % 2 == 1 && (j / square_size) % 2 == 1) )
            {
                checkerboard(i, j) = 0;
            }
            else
            {
                checkerboard(i, j) = 255;
            }
        }
    }

    normalizeExportImage(nullptr, "checkerboard", "png", checkerboard, 200, 200);

    return checkerboard;
}

MatrixXd addNoisetoImage(int width, int height, MatrixXd image_matrix)
{
    MatrixXd noise(height, width);
    noise.setRandom();                   // Fill the matrix with random values between -1 and 1
    noise = noise * 50.0;                // Scale the noise to be between -50 and 50
    image_matrix = image_matrix + noise; // Add the noise to the original image
    
    normalizeExportImage(nullptr, "noisy_image", "png", image_matrix, width, height);
    return image_matrix;
}

void computeLargestSingularEigenValues(MatrixXd singularValues, int numberOfValues){

    cout << "The " << numberOfValues << " largest singular values are: " << endl;
    for (int i = 0; i < numberOfValues; ++i) {
        cout << singularValues(i) << endl;
    }
}
 
double computeMeanSquareError(const MatrixXd &image1, const MatrixXd &image2)
{
    double mse = (image1 - image2).squaredNorm() / image1.size();
    return mse;
    
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
   
     // -- Task 4
    cout << "\n--------TASK 4---------- Not sure\n";
    task4_shiftIterativeSolver();
 
   // -- Task 5
    cout << "\n--------TASK 5----------\n";
    JacobiSVD<MatrixXd> svdA = SVD(image_matrix);

    // -- Task 6
    cout << "\n--------TASK 6----------\n";
    pair<MatrixXd, MatrixXd> cd40= createCD(svdA, 40);
    pair<MatrixXd, MatrixXd> cd80= createCD(svdA, 80);

      // -- Task 7
    cout << "\n--------TASK 7----------\n";
    MatrixXd compressedA40 = computeExportCompressedImage(cd40, "compressedA40");
    MatrixXd compressedA80 =computeExportCompressedImage(cd80, "compressedA80");

    // -- Task 8
    cout << "\n--------TASK 8----------\n";
    MatrixXd checkerboard=createCheckerBoard(200,25);

     // -- Task 9
    cout << "\n--------TASK 9----------\n";
    MatrixXd noisyCheckerBoard=addNoisetoImage(200, 200, checkerboard);

    // -- Task 10
    cout << "\n--------TASK 10----------\n";
    JacobiSVD<MatrixXd> svd_checkerboard = SVD(noisyCheckerBoard);
    VectorXd eigenvalues_checkerboard = svd_checkerboard.singularValues();
    computeLargestSingularEigenValues(eigenvalues_checkerboard, 2);

    // -- Task 11
    cout << "\n--------TASK 11----------\n";
    pair<MatrixXd, MatrixXd> cd5= createCD(svd_checkerboard, 5);
    pair<MatrixXd, MatrixXd> cd10= createCD(svd_checkerboard, 10);

    // -- Task 12
    cout << "\n--------TASK 12----------\n";
    MatrixXd compressedchb5 =computeExportCompressedImage(cd5, "compressedCheckerBoard5");
    MatrixXd compressedchb10 =computeExportCompressedImage(cd10, "compressedCheckerBoard10");

   // -- Task 13
   // Comparsing the compressed image with the original image
    cout << "\n--------TASK 13----------\n";
    double mse_compressed5 = computeMeanSquareError(checkerboard,compressedchb5);
    cout << "Mean Square Error between original checkerboard and compressed checkerboard with k=5: " << mse_compressed5 << endl;
    
    double mse_compressed10 = computeMeanSquareError(checkerboard,compressedchb10);
    cout << "Mean Square Error between original checkerboard and compressed checkerboard with k=10: " << mse_compressed10 << endl;

    double mse_noisy = computeMeanSquareError(checkerboard,noisyCheckerBoard);
    cout << "Mean Square Error between original checkerboard and noisy checkerboard: " << mse_noisy << endl;

    double mse_A40 = computeMeanSquareError(image_matrix,compressedA40);
    cout << "Mean Square Error between original image and compressed image with k=40: " << mse_A40 << endl;

     double mse_A80 = computeMeanSquareError(image_matrix,compressedA80);
    cout << "Mean Square Error between original image and compressed image with k=80: " << mse_A80 << endl;
    return 0;
}