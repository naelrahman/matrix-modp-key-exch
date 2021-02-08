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
    
    def __eq__(self, other):
        if (self.values == other.values):
            return True
        else:
            return False
        
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


# xk and xk_tuple square-and-multiply algorithm for efficiency.
def xk(self, k):
    k_binary = bin(k)
    value = MatrixModP(self.prime, *IDENTITY_MATRIX)
    
    # Starts at 2 to omit the "0b" part of binary representation.
    for digit in range(2, len(k_binary)):
        value *= value # Squares no matter what the binary digit is.
        if (k_binary[digit] == "1"): # Multiply by x when digit = 1.
            value *= self
    
    return value  

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
# (i) Generating M, H1, H2:
allCriteriaFit = False

while (allCriteriaFit == False):
    allCriteriaFit = False
    MATRIX_M_VALUES = []
    MATRIX_H1_VALUES = []
    MATRIX_H2_VALUES = []
    
    for x in range(9):
        MATRIX_M_VALUES.append(randbelow(PRIME_NUMBER))
        MATRIX_H1_VALUES.append(randbelow(PRIME_NUMBER))
        MATRIX_H2_VALUES.append(randbelow(PRIME_NUMBER))

    M = MatrixModP(PRIME_NUMBER, *MATRIX_M_VALUES)
    H1 = MatrixModP(PRIME_NUMBER, *MATRIX_H1_VALUES)
    H2 = MatrixModP(PRIME_NUMBER, *MATRIX_H2_VALUES)

    # If M is not invertible, redo generation:
    while (M.determinant() == 0):
        MATRIX_M_VALUES = [] # Clears matrix.
        for x in range(9):
            MATRIX_M_VALUES.append(randbelow(PRIME_NUMBER))

        M = MatrixModP(PRIME_NUMBER, *MATRIX_M_VALUES)
        print("Failed at M determinant")

    # If H is not invertible, redo generation:
    while (H1.determinant() == 0):
        MATRIX_H1_VALUES = []
        for x in range(9):
            MATRIX_H1_VALUES.append(randbelow(PRIME_NUMBER))
            
        H1 = MatrixModP(PRIME_NUMBER, *MATRIX_H1_VALUES)
        print("Failed at H1 determinant")

    while (H2.determinant() == 0):
        MATRIX_H2_VALUES = []
        for x in range(9):
            MATRIX_H2_VALUES.append(randbelow(PRIME_NUMBER))
            
        H2 = MatrixModP(PRIME_NUMBER, *MATRIX_H2_VALUES)
        print("Failed at H2 determinant")
    
    # Redoes generation of M and H1 or H2 if they are not invertible:
    MH = M * H1
    HM = H1 * M
    if MH.values == HM.values:
        print("Failed at H, M1 invertible")        
        continue 
    
    MH = M * H2
    HM = H2 * M
    if MH.values == HM.values:
        print("Failed at H, M2 invertible")                
        continue 
    
    # Ensures order of H^(p+1)(p^2+p+1) != 0
    orderCheckH1 = xk(H1, ((PRIME_NUMBER**2 + PRIME_NUMBER + 1) * (PRIME_NUMBER + 1)))
    orderCheckH2 = xk(H2, ((PRIME_NUMBER**2 + PRIME_NUMBER + 1) * (PRIME_NUMBER + 1)))
    
    if (orderCheckH1.values != IDENTITY_MATRIX and orderCheckH2.values != IDENTITY_MATRIX):
        print("Checks if identity matrix occurs") 
        allCriteriaFit = True

print("Finished generating M, H1, H2")

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
AliceCalculation = xk_tuple((M, (H1, H2)), m)
A = AliceCalculation[0]

print("Finished calculating A")
     
# STEP 3: Bob calculates and sends B.
BobCalculation = xk_tuple((M, (H1, H2)), n)
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
KeyCalculation = xk_tuple((M, (H1, H2)), m + n)
K = KeyCalculation[0]

print("Finished calculating K")


print("The prime number is set to", PRIME_NUMBER)
print("The chosen m, n values are", "m =", m, "n =", n)

print("\n------------")

print("\nGenerated M:")
M.output()

print("\nGenerated H1:")
H1.output()

print("\nGenerated H2:")
H2.output()

print("\n------------")

print("\nAlice's Key (K_A)")
K_A.output()

print("\nBob's Key (K_B)")
K_B.output()

print("\nShared Secret Key (K)")
K.output()

if (K_A == K):
    print("\nAlice's Key is equal to the Actual Key")
    
if (K_B == K):
    print("Bob's Key is equal to the Actual Key")

