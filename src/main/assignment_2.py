import numpy as np

# Question 1 || Using Nevill Method
# find 2nd degree, interpolating value for f(3.7)
#      x            f(x)
#     3.6          1.675
#     3.8          1.436
#     3.9          1.318

def nevilles_method(x_points, y_points, x):
  # must specify the matrix size (this is based on how many columns/rows you want)
  matrix = np.zeros((len(x_points), len(x_points)))
  # fill in value (just the y values because we already have x set)
  for counter, row in enumerate(matrix):
    row[0] = y_points[counter]
  # the end of the first loop are how many columns you have...
  num_of_points = len(x_points)

  # populate the matrix with the formula
  for i in range(1, num_of_points):
    for j in range(1, i+1):
      first_multiplication = (x- x_points[i-j]) * matrix[i][j-1]
      second_multiplication = (x - x_points[i]) * matrix[i-1][j-1]
      denominator = x_points[i] - x_points[i-j]
      coefficient = (first_multiplication - second_multiplication)/denominator
      matrix[i][j] = coefficient
    
  return matrix[num_of_points - 1][num_of_points - 1] #or coefficient

# Question 2 || Using Nevill Newton's Forward Method
# print out the polynomial approximations for degrees 1, 2, 3
#      x            f(x)
#     7.2          23.5492 
#     7.4          25.3913 
#     7.5          26.8224 
#     7.6          27.4589 
def divided_difference_table(x_points, y_points):
  size = len(x_points)
  matrix = np.zeros((size, size))
  result = []

  for i in range(size):
     matrix[i][0] = y_points[i]

  for i in range(1, size):
    for j in range(1, i + 1):
      matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (x_points[i] - x_points[i - j])

      if i == j:
        result.append(matrix[i][j])
  
  print(result)
  return matrix

# Question 3 || approximate f(7.3)
def get_approximate_result(matrix, x_points, value, start):
  reoccuring_x_span = 1
  reoccuring_px_result = start
    
  # we only need the diagonals...and that starts at the first row...
  for index in range(1, len(matrix)):
    polynomial_coefficient = matrix[index][index]

    # we use the previous index for x_points....
    reoccuring_x_span *= (value - x_points[index - 1])
        
    # get a_of_x * the x_span
    mult_operation = polynomial_coefficient * reoccuring_x_span

    # add the reoccuring px result
    reoccuring_px_result += mult_operation

  print(reoccuring_px_result)


# Question 4 || Hermite Polynomial Approximation Matrix
#      x            f(x)           f'(x)
#     3.6           1.675         -1.195
#     3.8           1.436         -1.188
#     3.9           1.318         -1.182 

np.set_printoptions(precision=7, suppress=True, linewidth=100)
def apply_div_dif(matrix: np.array):
  size = len(matrix)
  for i in range(2, size):
    for j in range(2, i+2):
      # skip if value is prefilled (we dont want to accidentally recalculate...)
      if j >= len(matrix[i]) or matrix[i][j] != 0:
        continue
            
      # get left cell entry
      left: float = matrix[i-1][j]
      # get diagonal left entry
      diagonal_left: float = matrix[i-1][j-1]
      # order of numerator is SPECIFIC.
      numerator: float = diagonal_left - left
      # denominator is current i's x_val minus the starting i's x_val....
      denominator = matrix[i][j] - matrix[i][i-j]
      # something save into matrix
      operation = numerator / denominator
      matrix[i][j] = operation
    
  return matrix

def hermite(x_points, y_points, derivative):
  # matrix size changes because of "doubling" up info for hermite 
  size = 2 * len(x_points)
  matrix = np.zeros((size, size + 1))

  # populate x values (make sure to fill every TWO rows)
  for i in range(size):
    matrix[i][0] = x_points[i // 2]
    matrix[i][1] = y_points[i // 2]

  # prepopulate y values (make sure to fill every TWO rows)
  for i in range(1, size, 2):
    matrix[i][2] = derivative[i // 2]

  # prepopulate with derivates (make sure to fill every TWO rows. starting row CHANGES.)
  for i in range(2, size):
    for j in range(2, i + 2):
      if matrix[i][j] != 0.:
        continue
            
      #difference formula
      matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (matrix[i][0] - matrix[i - j + 1][0])

  b = np.delete(matrix, -1, axis=1)
  print(b)

#Question 5
def cubic_spline_interpolation(x_points, y_points):
  size = len(x_points)
  matrix = np.zeros((size, size))
  matrix[0][0] = matrix[size - 1][size - 1] = 1

  for i in range(1, size - 1):
    index = i - 1
    for j in range(index, size, 2):
      matrix[i][j] = x_points[index + 1] - x_points[index]
      index += 1

  for i in range(1, size - 1):
    matrix[i][i] = 2 * ((x_points[i + 1] - x_points[i]) + (x_points[i] - x_points[i - 1]))

  print(np.matrix(matrix))
  print()

  spline_condition = np.zeros((size))

  for i in range (1, size - 1):
    first_term = (3 / (x_points[i + 1] - x_points[i])) * (y_points[i + 1] - y_points[i])
    second_term = (3 / (x_points[i] - x_points[i - 1])) * (y_points[i] - y_points[i - 1])
    spline_condition[i] = first_term - second_term

  print(np.array(spline_condition))
  print()

  print(np.array(np.linalg.solve(matrix, spline_condition)))

if __name__ == "__main__":
	# Question 1
  x_points = [3.6, 3.8, 3.9]
  y_points = [1.675, 1.436, 1.318]
  print(nevilles_method(x_points, y_points, 3.7))
  print()

  # Question 2
  # Task Two: print polynomial approximations for degrees 1, 2, 3
  x_points = [7.2, 7.4, 7.5, 7.6]
  y_points = [23.5492, 25.3913, 26.8224, 27.4589]
  divided_table = divided_difference_table(x_points, y_points)
  print()

  # Task Three: use the results from Task Two to approximate f(7.3)
  get_approximate_result(divided_table, x_points, 7.3, y_points[0])
  print()

  #Question 4
  hermite([3.6, 3.8, 3.9], [1.675, 1.436, 1.318], [-1.195, -1.188, -1.182])
  print()

  #Question 5
  cubic_spline_interpolation([2, 5, 8, 10], [3, 5, 7, 9])

  print()