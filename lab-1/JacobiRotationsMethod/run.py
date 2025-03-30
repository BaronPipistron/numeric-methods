import numpy as np


def is_symmetric_matrix(matrix) -> bool:
    matrix_size = matrix.shape[0]

    for i in range(matrix_size):
        for j in range(i + 1, matrix_size):
            if matrix[i, j] != matrix[j, i]:
                return False

    return True


def jacobi_rotations_method(matrix, epsilon: float):
    matrix_size = matrix.shape[0]

    matrix_cpy = np.copy(matrix)
    u_matrix = np.eye(matrix_size)

    iteration_counter = 0

    while True:
        iteration_counter += 1

        max_element_i, max_element_j = -1, -1
        max_element_abs_value = 0

        for i in range(matrix_size):
            for j in range(i + 1, matrix_size):
                if abs(matrix_cpy[i, j]) > max_element_abs_value:
                    max_element_i, max_element_j = i, j
                    max_element_abs_value = abs(matrix_cpy[i, j])

        phi = 0.5 * np.atan(
            (2 * matrix_cpy[max_element_i, max_element_j]) /
            (matrix_cpy[max_element_i, max_element_i] - matrix_cpy[max_element_j, max_element_j])
        )

        u_k_iteration = np.eye(matrix_size)
        u_k_iteration[max_element_i, max_element_i] = np.cos(phi)
        u_k_iteration[max_element_i, max_element_j] = -np.sin(phi)
        u_k_iteration[max_element_j, max_element_i] = np.sin(phi)
        u_k_iteration[max_element_j, max_element_j] = np.cos(phi)

        matrix_cpy = (u_k_iteration.T @ matrix_cpy) @ u_k_iteration

        u_matrix @= u_k_iteration

        t = 0
        for i in range(matrix_size):
            for j in range(i + 1, matrix_size):
                t += matrix_cpy[i, j] * matrix_cpy[i, j]

        t = np.sqrt(t)

        if t < epsilon:
            break

    eigenvalues = np.array([matrix_cpy[i, i] for i in range(matrix_size)])

    return eigenvalues, np.transpose(u_matrix), iteration_counter


def check_eigenvalues(matrix, eigenvalues, epsilon: float) -> bool:
    matrix_trace = sum(
        list(matrix[i][i] for i in range(matrix.shape[0]))
    )

    if abs(matrix_trace - sum(eigenvalues)) > epsilon:
        return False

    return True


def check_eigenvectors(
        matrix,
        eigenvalues,
        eigenvectors
) -> bool:
    for i in range(matrix.shape[0]):
        check_vector = (eigenvectors[i] @ matrix) / eigenvalues[i]

        if not np.allclose(check_vector, eigenvectors[i]):
            return False

    return True


if __name__ == '__main__':
    n = int(input('Input matrix size (one positive number): '))

    print('Input matrix:')
    input_matrix = np.empty((n, n), dtype=float)
    for row in range(n):
        input_matrix[row] = np.array(
            list(
                map(
                    float,
                    input().split()
                )
            )
        )

    eps = float(input('Input error rate: '))

    if not is_symmetric_matrix(input_matrix):
        print('ERROR: Jacobi Rotations Method not applicable - matrix is not symmetric!')
        exit(1)
    else:
        print('INFO: Jacobi Rotations Method applicable!')

    eigenvalues, eigenvectors, iteration_counter = jacobi_rotations_method(input_matrix, eps)

    print('Number of iterations: ', iteration_counter)
    print('\nEigenvalues:\n', eigenvalues)

    print('\nEigenvectors:\n', eigenvectors)
    print()

    if check_eigenvalues(input_matrix, eigenvalues, eps):
        print('Eigenvalues find correctly')
    else:
        print('Eigenvalues is not correct')

    if check_eigenvectors(input_matrix, eigenvalues, eigenvectors):
        print('Eigenvectors find correctly')
    else:
        print('Eigenvectors is not correct')

