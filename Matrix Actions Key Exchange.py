from math import gcd

from random import randint
from secrets import randbelow

# Note that the prime number should be a safe prime.
# The current prime is a 2000 bit number.
PRIME_NUMBER = int("""10045850546888500363341857765622433390255317048443
698327360730996384584773950711586086596475323993902797233883470790394194
018831434867898180891041375430671896508726694442987824141057899173376250
244281758576559881643143110828207143325627334593997352683778809319929255
772120459055406150435912157422236830704891980901048998096101770672922203
479101713092507042689334981405714581299534099154890607833310495144061448
203735644386469996712429901203439781034231264233355059817445403969916571
063605224058329470399818911447991765712527069708623420044248954447465956
0583354052797579309573507121265302226528942789519""".replace("\n", ""))

IDENTITY_MATRIX = [1, 0, 0, 0, 1, 0, 0, 0, 1]
ZERO_MATRIX = [0, 0, 0, 0, 0, 0, 0, 0, 0]

# MatrixModP: 3x3 matrix with entries in the range [0, PRIME_NUMBER)
class MatrixModP():
    def __init__(self, prime, *numbers):
        self.values = ZERO_MATRIX[:]
        self.prime = prime
        
        counter = 0
        for x in numbers: # Adjust values to user's inputs.
            self.values[counter] = x % self.prime
            counter += 1
            
    def __add__(self, other):
        sums_values = []
        for x in range(9):
            toInsert = self.values[x] + other.values[x]
            toInsert = toInsert % self.prime
            sums_values.append(toInsert)
        
        return MatrixModP(self.prime, *sums_values)
    
    def __mul__(self, other):
        product_values = []
        
        for x in range(3): # Loops through first matrix' rows.
            for y in range(3): # Loops through second matrix' columns.
                toInsert = 0
                for index in range(3):
                    toInsert += (self.values[(x * 3) + index] *
                        other.values[y + (index * 3)])
                    # self.values part: "x * 3" determines what row to
                    # look at, index loops through that row.
                    
                    # other.values part: "y" is the first row value,
                    # index * 3 cycles through the row
                
                toInsert = toInsert % self.prime    
                product_values.append(toInsert)
        
        return MatrixModP(self.prime, *product_values)
    
    def __pow__(self, k):
        k_binary = bin(k)
        value = MatrixModP(self.prime, *IDENTITY_MATRIX)
        
        # Starts at 2 to omit the "0b" part of binary representation.
        for digit in range(2, len(k_binary)):
            value *= value # Squares no matter what the binary digit is.
            if (k_binary[digit] == "1"): # Multiply by x when digit = 1.
                value *= self
        
        return value
    
    def __eq__(self, other):
        if (self.values == other.values):
            return True
        else:
            return False
    
    def change_val(self, index, new):
        self.prime[index] = new
        
    def output(self):
        start = 0
        end = 3
        
        for x in range(3):
            print(self.values[start:end])
            start = end
            end += 3
            
    def determinant(self):
        calculation = (self.values[0]*self.values[4]*self.values[8]
            + self.values[1]*self.values[5]*self.values[6]  
            + self.values[2]*self.values[3]*self.values[7]  
            - self.values[0]*self.values[5]*self.values[7]  
            - self.values[1]*self.values[3]*self.values[8]  
            - self.values[2]*self.values[4]*self.values[6]) 
        
        adjusted = (calculation % self.prime)
        return adjusted

# Procedure for creating matrix M:
def GenerateM():
    MATRIX_M_VALUES = []
    for x in range(9):
        MATRIX_M_VALUES.append(randbelow(PRIME_NUMBER))
    return MatrixModP(PRIME_NUMBER, *MATRIX_M_VALUES)

# Procedure for creating matrices H_1 and H_2.
def GenerateH():
    DIAGONAL_MATRIX = [0, 0, 0, 0, -1, 0, 0, 0, -1]
    for x in [4, 8]:
        while (DIAGONAL_MATRIX[x] == -1):
            toInsert = randbelow(PRIME_NUMBER)
            if (toInsert**2 % PRIME_NUMBER != 1):
                DIAGONAL_MATRIX[x] = toInsert
    
    D = MatrixModP(PRIME_NUMBER, *DIAGONAL_MATRIX)
    
    S_VALUES = [randbelow(PRIME_NUMBER) for i in range(9)]
    S = MatrixModP(PRIME_NUMBER, *S_VALUES)
    
    while (S.determinant == 0):
        S.change_val(0, randbelow(PRIME_NUMBER))
    
    return InverseMatrix(S) * D * S
    

def InverseMatrix(A):
    (L, U) = LUDecomp(A)
    
    INV_MAT_COLUMNS = ([1, 0, 0], [0, 1, 0], [0, 0, 1])
    col_1 = LUSolve(L, U, INV_MAT_COLUMNS[0])
    col_2 = LUSolve(L, U, INV_MAT_COLUMNS[1])
    col_3 = LUSolve(L, U, INV_MAT_COLUMNS[2])
    
    M_INV_VAL = [col_1[0], col_2[0], col_3[0],
        col_1[1], col_2[1], col_3[1],
        col_1[2], col_2[2], col_3[2]]

    return MatrixModP(PRIME_NUMBER, *M_INV_VAL)
    

def LUDecomp(A):
    A_VALS = A.values[:]
    U_VALS = [A_VALS[0], -1, -1, 0, A_VALS[4], -1, 0, 0, A_VALS[8]]
    L_VALS = [1, 0, 0, -1, 1, 0, -1, -1, 1]
    
    for k in range(3):
        U_VALS[(k*3)+k] = A_VALS[(k*3)+k]
        for i in range(k+1, 3):
            L_VALS[(i*3)+k] = (A_VALS[(i*3)+k] * pow(U_VALS[(k*3)+k], -1, PRIME_NUMBER))
            U_VALS[(k*3)+i] = A_VALS[(k*3)+i]
        for i in range(k+1, 3):
            for j in range(k+1, 3):
                A_VALS[(i*3)+j] -= L_VALS[(i*3)+k] * U_VALS[(k*3)+j]

    return (MatrixModP(PRIME_NUMBER, *L_VALS), MatrixModP(PRIME_NUMBER, *U_VALS))

def LUSolve(L, U, b):
    y = [1, 1, 1]
    x = [1, 1, 1]
    for i in range(3):
        toSub = 0
        for j in range(0, i):
            toSub += (L.values[(i*3)+j]*y[j])
        y[i] = b[i] - toSub
    
    for i in range(2, -1, -1):
        toSub = 0
        for j in range(i+1, 3):
            toSub += (U.values[(i*3)+j]*x[j])
        x[i] = (y[i] - toSub) * pow(U.values[(i*3)+i], -1, PRIME_NUMBER)
    
    return x

# Precondition: "self" should be a 2-tuple with a matrix and tuple.
#   At the beginning of the procedure, this is ((M, (H1, H2), k)
def xk_tuple(self, k):
    k_binary = bin(k)
    value = self
    
    # The reason this for-loop starts with 3 instead of 2 is because
    # value is initalized as self, which erases the need for the first
    # multiplication of the square-and-multiply method.
    for digit in range(3, len(k_binary)): 
        value = semidirect_product(value, value)[:]
        if (k_binary[digit] == "1"):
            value = semidirect_product(value, self)[:]
    
    return value  


# Semidirect product, as described in Section 2.
#   How firstTuple and secondTuple are broken down:
#       The second element of each of these tuples is another 2-tuple,
#       that contains (H1, H2) to some power..
def semidirect_product(firstTuple, secondTuple):
    resultingFirstMatrix = ((secondTuple[1][0] * firstTuple[0] * secondTuple[1][1])
        + secondTuple[0])
    
    resultingH1 = (firstTuple[1][0] * secondTuple[1][0])
    resultingH2 = (firstTuple[1][1] * secondTuple[1][1])
    
    return [resultingFirstMatrix, (resultingH1, resultingH2)]



# ----------------------------------------------------------------------
#                       PROTOCOL DESCRIPTION
# ----------------------------------------------------------------------


# STEP 1:
# (i) Generating M, H_1, H_2:
M = GenerateM()
H_1 = GenerateH()
H_2 = GenerateH()

for H in (H_1, H_2):
    while (M * H == H * M):
        M = GenerateM()
        print("uh oh stinky")
        
print("Finished generating M, H_1, H_2")

# (ii) Alice and Bob's private selections.
#   Range of m and n, based on the magnitude of q.
MAGNITUDE = len(str(PRIME_NUMBER))
if str(PRIME_NUMBER)[0] == "1": # Applies when mag of p != mag of q
    MAGNITUDE -= 1

LOWEST_EXPONENT = 10**(MAGNITUDE - 1)
HIGHEST_EXPONENT = 10**(MAGNITUDE)


m = randint(LOWEST_EXPONENT, HIGHEST_EXPONENT)
n = randint(LOWEST_EXPONENT, HIGHEST_EXPONENT)

#   Redoes generation of m, n until both are relatively prime to p - 1.
while (gcd(m, PRIME_NUMBER - 1) != 1 or gcd(n, PRIME_NUMBER - 1) != 1):
    m = randint(LOWEST_EXPONENT, HIGHEST_EXPONENT)
    n = randint(LOWEST_EXPONENT, HIGHEST_EXPONENT)

print("Finished generating m and n")

# STEP 2: Alice calculates and sends A.
AliceCalculation = xk_tuple((M, (H_1, H_2)), m)
A = AliceCalculation[0]

print("Finished calculating A")
     
# STEP 3: Bob calculates and sends B.
BobCalculation = xk_tuple((M, (H_1, H_2)), n)
B = BobCalculation[0]

print("Finished calculating B")

# STEP 4: Alice retrieves K_A by computing (B, x)* (A, (H1^m, H2^m)).
AliceKey = semidirect_product(BobCalculation, AliceCalculation)
K_A = AliceKey[0]

print("Finished calculating K_A")

# STEP 5: Bob retrieves K_B by computing (A, y)* (B, (H1^n, H2^n)).
BobKey = semidirect_product(AliceCalculation, BobCalculation)
K_B = BobKey[0]

print("Finished calculating K_B")

# STEP 6: Calculation of secret key K to confirm key generation works.
KeyCalculation = xk_tuple((M, (H_1, H_2)), m + n)
K = KeyCalculation[0]

print("Finished calculating K")


print("The prime number is set to", PRIME_NUMBER)
print("The chosen m, n values are", "m =", m, "n =", n)

print("\n------------")

print("\nGenerated M:")
M.output()

print("\nGenerated H1:")
H_1.output()

print("\nGenerated H2:")
H_2.output()

print("\n------------")

print("\nAlice's Key (K_A)")
K_A.output()

print("\nBob's Key (K_B)")
K_B.output()

print("\nShared Secret Key (K)")
K.output()

# Checks to see if the keys are equal to the actual keys.
if (K_A == K):
    print("\nAlice's Key is equal to the Actual Key")
    
if (K_B == K):
    print("Bob's Key is equal to the Actual Key")
