EPS = 1e-12

#Chuyen tat ca phan tu xap xi 0 thanh 0
#tranh sai sot trong tinh rank
def is_zero_row(row):
    return all(abs(x) < EPS for x in row)

def rank_and_basis(A):
    n = len(A)
    m = len(A[0])
    matrix = [row[:] for row in A]

    row = 0
    pivot_cols = []
    for col in range(m):
        pivot = None
        for r in range(row, n):
            if abs(matrix[r][col]) > EPS:
                pivot = r
                break

        if pivot is None:
            continue

        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        pivot_cols.append(col)
        for r in range(row + 1, n):
            factor = matrix[r][col] / matrix[row][col]
            for c in range(col, m):
                matrix[r][c] -= factor * matrix[row][c]

        row += 1

    rank = 0
    for r in matrix:
        if not is_zero_row(r):
            rank += 1

    #Lay khong gian tu nhung cot khac 0
    cols_space = []
    for col in pivot_cols:
        colum = []
        for i in range(n):
            colum.append(A[i][col])
        cols_space.append(colum)

    #dong khac 0
    rows_space = []
    for r in matrix:
        if not is_zero_row(r):
            rows_space.append(r)
    
    #Lay cac vi tri dong co bien tu do
    free_vars = []
    for j in range(m):
        found = False
        for pc in pivot_cols:
            if j == pc:
                found = True
                break
        if not found:
            free_vars.append(j)


    null_sppace = []
    for free in free_vars:
        #cho 1 bien tu do ban 1 con lai bang 0
        x = [0.0] * m
        x[free] = 1.0

        for i in range(rank - 1, -1, -1):
            pivot_col = pivot_cols[i]

            sum_ax = 0.0
            for j in range(pivot_col + 1, m):
                sum_ax += matrix[i][j] * x[j]

            if abs(matrix[i][pivot_col]) < EPS:
                x[pivot_col] = 0
            else:
                x[pivot_col] = -sum_ax / matrix[i][pivot_col]

        null_sppace.append(x)

    return rank, rows_space, cols_space, null_sppace