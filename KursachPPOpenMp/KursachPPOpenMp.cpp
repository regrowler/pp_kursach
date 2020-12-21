#include <iostream>
#include <vector>
#include <iostream>
#include "omp.h"
#include "linear_system_input.h"
#include <ctime>/*
#pragma comment(lib,"C:\\Users\\AndreiKoloskov\\Downloads\\mpfr_mpir_x86_x64_msvc2010\\mpir\\dll\\x64\\Debug\\mpir.lib")
#include "C:\\Users\\AndreiKoloskov\\Downloads\\mpfr_mpir_x86_x64_msvc2010\\mpir\\dll\\x64\\Debug\\mpir.h"
#include "C:\\Users\\AndreiKoloskov\\Downloads\\mpfr_mpir_x86_x64_msvc2010\\mpir\\dll\\x64\\Debug\\gmp.h"*/
using namespace std;

vector<vector<double>> copy_matrix(vector<vector<double>> arr, int N) {
    vector<vector<double>> u = {};
    for (int i = 0; i < N; i++) {
        u.push_back({});
        for (int j = 0; j < N; j++) {
            u[i].push_back(arr[i][j]);
        }
    }
    return u;
}


double determinant1(std::vector<std::vector<double>>& matrix) {
    int N = static_cast<int>(matrix.size());
    double det = 1;

    for (int i = 0; i < N; ++i) {

        double pivotElement = matrix[i][i];
        int pivotRow = i;
        for (int row = i + 1; row < N; ++row) {
            if (std::abs(matrix[row][i]) > std::abs(pivotElement)) {
                pivotElement = matrix[row][i];
                pivotRow = row;
            }
        }
        if (pivotElement == 0.0) {
            return 0.0;
        }
        if (pivotRow != i) {
            matrix[i].swap(matrix[pivotRow]);
            det *= -1.0;
        }
        det *= pivotElement;

        for (int row = i + 1; row < N; ++row) {
            for (int col = i + 1; col < N; ++col) {
                matrix[row][col] -= matrix[row][i] * matrix[i][col] / pivotElement;
            }
        }
    }
    std::cout << "determinant gauss" << det << std::endl;
    return det;
}
//
//__mpq_struct determinant3(std::vector<std::vector<double>>& matrix) {
//    int N = static_cast<int>(matrix.size());
//    mpq_t det;
//    mpq_init(det);
//    __gmpq_set_d(det, 0.0);
//
//    for (int i = 0; i < N; ++i) {
//
//        double pivotElement = matrix[i][i];
//        int pivotRow = i;
//        for (int row = i + 1; row < N; ++row) {
//            if (std::abs(matrix[row][i]) > std::abs(pivotElement)) {
//                pivotElement = matrix[row][i];
//                pivotRow = row;
//            }
//        }
//        if (pivotElement == 0.0) {
//            __gmpq_set_d(det, 0.0);
//            return det[0];
//        }
//        if (pivotRow != i) {
//            matrix[i].swap(matrix[pivotRow]);
//            mpq_t temp;
//            mpq_init(temp);
//            __gmpq_set_d(temp, 0.0);
//            mpq_mul(det, det, temp);
//            return temp[0];
//        }
//        mpq_t temp;
//        mpq_init(temp);
//        __gmpq_set_d(temp, pivotElement);
//        mpq_mul(det, det, temp);
//
//        for (int row = i + 1; row < N; ++row) {
//            for (int col = i + 1; col < N; ++col) {
//                matrix[row][col] -= matrix[row][i] * matrix[i][col] / pivotElement;
//            }
//        }
//    }
//    std::cout << "determinant gauss" << det << std::endl;
//    return det[0];
//}

void swap_column(vector<vector<double>>& matrix,vector<double> newColumn, size_t i, size_t j) {
    for (size_t t = 0; t < j; t++)
    {
        matrix[t][i] = newColumn[t];
    }
    int tt = 0;
}

std::vector<double> cramer_solving_omp(std::vector<std::vector<double> >& matrix,std::vector<double>& free_term_column, int n) {
    vector<vector<double>> copiedMatrix = copy_matrix(matrix, n);
    double mainDet = determinant1(copiedMatrix);
    if (std::abs(mainDet) < 0.0001)
        cout << "no desicion" << endl;
    std::vector<double> solution(n);
        
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        std::vector<std::vector<double> > private_matrix = copy_matrix(matrix,n);
        swap_column(private_matrix, free_term_column, i, n);
        solution[i] = determinant1(private_matrix) / mainDet;
    }
    return solution;
}

//std::vector<__mpq_struct> cramer_solving_omp_mpir(std::vector<std::vector<double> >& matrix, std::vector<double>& free_term_column, int n) {
//    vector<vector<double>> copiedMatrix = copy_matrix(matrix, n);
//    double mainDet = determinant1(copiedMatrix);
//    if (std::abs(mainDet) < 0.0001)
//        cout << "no desicion" << endl;
//    vector<__mpq_struct> solution(n);
//    #pragma omp parallel for
//    for (int i = 0; i < n; ++i) {
//        std::vector<std::vector<double> > private_matrix = copy_matrix(matrix, n);
//        swap_column(private_matrix, free_term_column, i, n);
//        mpq_t asd = { determinant3(private_matrix) };
//        solution[0] = asd[0];
//        //mpq_add(solution[0],iqw,iqw)
//        //solution[i] = determinant1(private_matrix) / mainDet;
//    }
//    return solution;
//}

//std::vector<__mpq_struct> cramer_solving_mpir(std::vector<std::vector<double> >& matrix, std::vector<double>& free_term_column, int n) {
//    vector<vector<double>> copiedMatrix = copy_matrix(matrix, n);
//    double mainDet = determinant1(copiedMatrix);
//    if (std::abs(mainDet) < 0.0001)
//        cout << "no desicion" << endl;
//    vector<__mpq_struct> solution(n);
//    for (int i = 0; i < n; ++i) {
//        std::vector<std::vector<double> > private_matrix = copy_matrix(matrix, n);
//        swap_column(private_matrix, free_term_column, i, n);
//        mpq_t asd = { determinant3(private_matrix) };
//        solution[0] = asd[0];
//        //mpq_add(solution[0],iqw,iqw)
//        //solution[i] = determinant1(private_matrix) / mainDet;
//    }
//    return solution;
//}

std::vector<double> cramer_solving(std::vector<std::vector<double> >& matrix, std::vector<double>& free_term_column, int n) {
    vector<vector<double>> copiedMatrix = copy_matrix(matrix, n);
    double mainDet = determinant1(copiedMatrix);
    if (std::abs(mainDet) < 0.0001)
        cout << "no desicion" << endl;
    std::vector<double> solution(n);
    for (int i = 0; i < n; ++i) {
        std::vector<std::vector<double> > private_matrix = copy_matrix(matrix, n);
        swap_column(private_matrix, free_term_column, i, n);
        solution[i] = determinant1(private_matrix) / mainDet;
    }
    return solution;
}
int main() {
    omp_set_num_threads(8);
    srand(time(0));
    LinearSystemInput input;
    input.read_equations_from_file("C:\\Users\\AndreiKoloskov\\Documents\\kramer\\input.txt");
    int matrix_rank = input.get_rank();
    int buffer_size = matrix_rank * matrix_rank;
    vector<vector<double>> matrix = input.storage;
    vector<double> freTerms = {};
    for (size_t i = 0; i < matrix_rank; i++)
    {
        for (size_t j = 0; j <= matrix_rank; j++)
        {
            //cout << matrix[i][j] << " ";
            if (j== matrix_rank)
            {
                freTerms.push_back(matrix[i][j]);
                matrix[i].pop_back();
            }
        }
        //cout << endl;
    }
    vector<vector<double>>m1 = copy_matrix(matrix, matrix_rank);
    vector<vector<double>> m2 = copy_matrix(matrix, matrix_rank);
    
    cout << "solving with sequintial" << endl;
    double start1 = (double)clock() / 1000;
    std::vector<double> solution1 = cramer_solving(m1, freTerms, matrix_rank);
    double end1 = (double)clock() / 1000;
    
    cout << "sequintial method result" << endl;
    for (size_t i = 0; i < matrix_rank; i++)
    {
        //gmp_printf("the rational is: %Qd\n", a);
        std::cout << "X" << i + 1 << " = " << solution1[i] << std::endl;
    }
    cout << "time of solving a system of equations with dimension " << matrix_rank << " with sequintial: " << end1 - start1 << " seconds" << endl;


    cout << "solving with omp" << endl;
    double start2 = (double)clock() / 1000;
    std::vector<double> solution2 = cramer_solving_omp(m1, freTerms, matrix_rank);
    double end2 = (double)clock() / 1000;
    cout << "omp method result" << endl;
    for (size_t i = 0; i < matrix_rank; i++)
    {
        std::cout << "X" << i + 1 << " = " << solution2[i] << std::endl;
    }
    cout << "time of solving a system of equations with dimension " << matrix_rank << " with openMp: " << end2 - start2 << " seconds" << endl;
   
  
}

