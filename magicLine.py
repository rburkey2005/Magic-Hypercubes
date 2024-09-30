#!/usr/bin/env python3
program = "magicLine.py"
revision = 0.14
"""
Copyright:  None.  This software is declared to be in the Public Domain by its 
            author, Ronald S. Burkey.
Purpose:    This is a replacement for magicHypercube.py.
            The idea is to construct a magic WxWxW cube from a magic WxW 
            square, a magic WxWxWxW hypercube from a magic WxWxW cube, and
            so on, using an algorithm I've invented.
History:    2022-07-12 RSB  Began
            2022-07-14 RSB  Added W=14, 18, 22.  Added some more ways of 
                            constructing p and q that filled in more gaps as
                            well.
            2022-07-15 RSB  Tries to track free memory and abort if there's
                            not enough.  Should allow running with very large
                            maxN without danger of crashing system. 
            2022-07-19 RSB  molsPower2() improved to produce self-complementary
                            MOLS P and Q, rather than just diagonal MOLS.
            2022-07-27 RSB  Removed W=10 from among the hard-coded cases to
                            an algorithmically-generated case, as part of the 
                            introduction of algorithmic generation of MOLS
                            for the case W=q+r, where q, r are prime and 
                            q = 2r + 1.
            2022-07-28 RSB  Realized that r did not have to be prime in the
                            scenario introduced yesterday, so enabled that
                            and cleaned up the implementation somewhat.  The
                            --debug messages about whether P and Q are
                            diagonal latin squares were wrong and have been
                            corrected.  Added --test-mols.
            2022-07-29 RSB  Added --qr-test
            2022-07-30 RSB  (0.11) Added --eng-test.  
            2022-07-31 RSB  Relaxed conditions on q+r construction.
            2022-08-01 RSB  (0.12) I hadn't implemented the change from 
                            yesterday quite correctly.  But a correct
                            implementation doesn't seem to change the results.
                            (0.13) Replaced the hardcoded p14 and q14 by a
                            hardcoded self-orthogonal transversal p14.
            2022-08-03 RSB  (0.14) Usage of hardcoded MOLS P and Q cleaned up.
                            Hardcoded MOLS moved into a separate imported
                            file, hardcodedMOLS.py.
"""

import sys
import copy
import os
import signal
import time
# The file hardcodedMOLS.py should be in the same folder as magicLine.py.
from hardcodedMOLS import hardcodedMOLS

# Return max memory usage (VmPeak), in a float GB.  The function is designed
# for Linux, and may return 0 elsewhere.
_FIELDS = ['VmRSS', 'VmHWM', 'VmSize', 'VmPeak']
def getVmpeakGB():
    try:
        # read in process info
        with open('/proc/self/status', 'r') as file:
            lines = file.read().split('\n')
        # check all process info fields
        for line in lines:
            if ':' in line:
                name, val = line.split(':')
    
                # collect relevant memory fields
                if name == "VmPeak":
                    value = int(val.strip().split(' ')[0])  # strip "kB"
                    return value / 1000000.0  # convert to GB
    except:
        pass
    return 0.0

startTime = time.time()
# Print a "report" message.
def printReport(msg, D):
    print("Report W=%d D=%d: %s (Time=%0.2f sec, Mem=%0.2f GB)" \
            % (W, D, msg, time.time()-startTime, getVmpeakGB()))

print()
print("Program", program, "revision", revision)

# Set defaults for command-line settings.
maxN = 3    # Maximum dimension N for which we want to generate a magic N-cube.
            # Must be > 2.
W = 3       # Width of the hypercubes.  I.e., the magic squares are WxW, the
            # magic cubes are WxWxW, and so on.  Must be > 2.
debug = False
onlyPartial = False
quiet = False
python = False
maxMemGB = 4.0
testMOLS = False
qrTest = 0
engTest = ""
try:
    error = False
    for param in sys.argv[1:]:
        fields = param.split("=")
        switch = fields[0]
        value = ""
        if len(fields) > 1:
            value = fields[1]
        if switch == "--dimension":
            maxN = int(value)
        elif switch == "--width":
            W = int(value)
        elif switch == "--debug":
            debug = True
        elif switch == "--partial":
            onlyPartial = True
        elif switch == "--quiet":
            quiet = True
        elif switch == "--python":
            python = True
        elif switch == "--max-memory":
            maxMemGB = float(value)
        elif switch == "--test-mols":
            testMOLS = True
        elif switch == "--qr-test":
            qrTest = int(value)
        elif switch == "--eng-test":
            engTest = value
        else:
            print("=== Unknown switch:", param, "===")
            error = True
except:
    error = True
if W < 3 or maxN < 2:
    error = True
if error:
    print("This program prints N-dimensional magic WxWx...xW hypercubes for")
    print("N=3,...,maxN. Usage:")
    print("     %s [OPTIONS]" % program)
    print("The available OPTIONS are:")
    print("--dimension=N     Sets maxN (default 3, mininum 2).")
    print("--width=N         Sets W (default 3, minimum 3). Not all values")
    print("                  of W are supported.  W=6 is impossible.")
    print("--partial         Check for partially- rather than fully-magic.")
    print("--quiet           Do not print the generated cubes/hypercubes.")
    print("                  This is useful for high dimensions or high W,")
    print("                  when the results of the checks are useful")
    print("                  rather than the result of the construction.")
    print("--python          Outputs result in the form of a Python list of")
    print("                  lists of ..., rather than human readable.")
    print("--max-memory=F    The program will try to exit gracefully if")
    print("                  its memory usage exceeds F GB.  The default is")
    print("                  4.0.")
    print("--debug           Enables extra debugging tests and messages.")
    print("--test-mols       Simply generates and tests MOLS, nothing more.")
    print("                  Ignores all other switches, including --width.")
    print("--qr-test=R       Enables a test mode for W=q+r constructions.")
    print("--eng-test=L      Adds parameters for whatever engineering tests")
    print("                  I feel like permforming at any given time.")
    print()
    sys.exit(1)

print("Width", W, "======================================================")

if W == 6:
    printReport("Impossible", 2)
    sys.exit(1)

smallPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
               59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
               127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
               191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 
               257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 
               331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 
               401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
               467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 
               563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 
               631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 
               709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 
               797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 
               877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 
               967, 971, 977, 983, 991, 997]  

# Starting values (1-dimensional) of the "magic" and "adder" objects.  These
# are the objects that will be converted into magic/adder squares, then cubes,
# then hypercubes.  Thus m1 is a "magic line" and a1 is a "latin line"; for
# both of which, [0, 1, ..., W-1] serves fine.
m = list(range(W))
a = list(range(W))

# Check memory use and abort if excessive.
def checkExcessiveMemory(W, N):
    vmPeak = getVmpeakGB()
    if vmPeak > maxMemGB:
        printReport("Excessive memory usage %0.2f > %0.2f GB" \
                        % (vmPeak, maxMemGB), N)
        sys.stdout.flush()
        time.sleep(1)
        sys.exit(1)

# Predict expected memory use and abort if excessive.
def predictExcessiveMemory(W, N):
    vmPeak = getVmpeakGB()
    predicted = 1.1 * W * vmPeak
    if predicted > maxMemGB:
        printReport("Predicted excessive memory usage %0.2f > %0.2f GB" \
                        % (predicted, maxMemGB), N)
        sys.stdout.flush()
        time.sleep(1)
        sys.exit(1)

# Transpose a square matrix, creating a new one.
def transposeMatrix(m):
    t = []
    for row in range(len(m)):
        trow = []
        for col in range(len(m)):
            trow.append(m[col][row])
        t.append(trow)
    return t

# Returns MOLS P and Q for W=2**(2*K+1). 
def molsPower2(W):
    twoTo2K = W // 2
    twoTo2Kminus1 = twoTo2K - 1
    r = []
    for i in range(W):
        r.append(i ^ ((i // twoTo2K) * twoTo2Kminus1))
    c = []
    for j in range(W):
        c.append((j // 2) + twoTo2K * ((j % 2) ^ (j // twoTo2K)))
    P = []
    for i in range(W):
        P.append([0]*W)
    for i in range(W):
        for j in range(W):
            P[i][j] = r[i] ^ c[j]
    return P, transposeMatrix(P)

# Returns MOLS P and Q for W an odd number not divisible by 3.
def oddNotDivisibleBy3(W):
    P = []
    for row in range(W):
        P.append([])
        for col in range(W):
            P[row].append((2*col - row) % W)
    return P, transposeMatrix(P)

# Generate a latin square of order W=p+r from two latin squares of order
# p and r, where p = 2r+1 is prime. The input 'reverse' should be False when
# constructing P from Pp and Pr, but should be True when constructing
# Q from Qp and Qr.
def oddPlusOdd(Pp, Pr, reverse):
    p = len(Pp)
    r = len(Pr)
    if qrTest == 3:
        qPrime = p // 3
        qPrimePrime = p % 3
    #if p != 2 * r + 1:
    #    return []
    W = p + r
    P = []
    for i in range(W): 
        P.append([0]*W)
    for i in range(p):
        for j in range(p):
            P[i][j] = Pp[i][j]
    for i in range(p):
        for k in range(r):
            if qrTest == 3 and engTest == "self":
                if k == 0:
                    sub = P[i][(p - 1 + i) % p]
                elif k == 1:
                    sub = P[i][(2*qPrime+qPrimePrime+i) % p]
                elif k == 2:
                    sub = P[i][(qPrime+1+i) % p]
            elif reverse:
                sub = P[i][(i + 1 + k) % p]
            else:
                sub = P[i][(p - 1 + i - k) % p]
            if engTest == "sub8" and W == 14 and i == 0 and k == 0:
                sub = 8
            P[i][p + k] = sub
            index = P[i][:p].index(sub)
            P[p + k][index] = P[i][index]
    for k in range(r):
        for l in range(r):
            P[p + k][p + l] = p + Pr[k][l]
    for i in range(p):
        for k in range(r):
            index = P[i][:p].index(P[i][p + k])
            P[i][index] = p + k
    return P
    
# Get a pair of mutually-orthogonal latin squares (MOLS) p,q for the chosen 
# value of W. If possible, they should have diagonals and antidiagonals of 
# all distinct values; failing that, the diagonals and antidiagonals should
# sum to W(W-1)/2; failing that ... well, I'll take whatever I can get.
# Return [],[] upon failure.  In some cases I don't have a satisfactory
# way of generating them yet, so they're hard-coded instead. In the case of
# W=6, it's simply theoretically impossible.  In still other cases, I
# have no method available as of yet.
def getPermutator(W):

    if W < 3:
        return [], []

    # Factor W = w * 2**pow2.
    pow2 = 0
    w = W
    while w%2 == 0:
        pow2 += 1
        w = w//2
    
    # This function takes a divide-and-conquer approach.  If W can be 
    # factored as the product of two factors, for each of which good
    # MOLS are known, then the MOLS for W can be constructed
    # from the MOLS of the factors.  FYI: Sometimes getForDivisor(W//divisor)
    # works better than getForDivisor(divisor); it's not totally symmetric,
    # in the sense that you the MOLS of the divisors to be self-complementary,
    # but it's easy to make a mistake and choose a divisor that doesn't
    # have self-complementary MOLS such as (at this writing) 10, 14, 18, or 22.
    def getForDivisor(divisor):
        pD, qD = getPermutator(divisor)
        if pD == []:
            return pD, qD
        dW = W // divisor
        pSmall,qSmall = getPermutator(dW)
        if pSmall == []:
            return pSmall, qSmall
        p = []
        q = []
        for rowBig in range(divisor):
            for rowSmall in range(dW):
                row = len(p)
                p.append([0]*W)
                q.append([0]*W)
                for colBig in range(divisor):
                    for colSmall in range(dW):
                        col = colBig * dW + colSmall
                        p[row][col] = pSmall[rowSmall][colSmall] + \
                                dW * pD[rowBig][colBig]
                        q[row][col] = qSmall[rowSmall][colSmall] + \
                                dW * qD[rowBig][colBig]
        return p, q
    
    # Some engineering-test cases
    if qrTest > 0 and W > 2 * qrTest \
            and W % 2 == 0: # and W // 2 in smallPrimes:
        small = qrTest
        big = W - small
        if big != 15 or engTest != "natalia":
            Pbig, Qbig = getPermutator(big)
        else:
            print("Here!")
            # A pair of transversal MOLS for W=15, Natalia Makarova.
            Pbig = [[0, 6, 11, 9, 8, 2, 14, 5, 4, 12, 3, 10, 13, 1, 7],
                    [2, 1, 7, 12, 10, 9, 3, 14, 6, 5, 13, 4, 11, 0, 8],
                    [1, 3, 2, 8, 13, 11, 10, 4, 14, 7, 6, 0, 5, 12, 9],
                    [13, 2, 4, 3, 9, 0, 12, 11, 5, 14, 8, 7, 1, 6, 10],
                    [7, 0, 3, 5, 4, 10, 1, 13, 12, 6, 14, 9, 8, 2, 11],
                    [3, 8, 1, 4, 6, 5, 11, 2, 0, 13, 7, 14, 10, 9, 12],
                    [10, 4, 9, 2, 5, 7, 6, 12, 3, 1, 0, 8, 14, 11, 13],
                    [12, 11, 5, 10, 3, 6, 8, 7, 13, 4, 2, 1, 9, 14, 0],
                    [14, 13, 12, 6, 11, 4, 7, 9, 8, 0, 5, 3, 2, 10, 1],
                    [11, 14, 0, 13, 7, 12, 5, 8, 10, 9, 1, 6, 4, 3, 2],
                    [4, 12, 14, 1, 0, 8, 13, 6, 9, 11, 10, 2, 7, 5, 3],
                    [6, 5, 13, 14, 2, 1, 9, 0, 7, 10, 12, 11, 3, 8, 4],
                    [9, 7, 6, 0, 14, 3, 2, 10, 1, 8, 11, 13, 12, 4, 5],
                    [5, 10, 8, 7, 1, 14, 4, 3, 11, 2, 9, 12, 0, 13, 6],
                    [8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 14]]
            Qbig = [[0, 10, 7, 2, 5, 11, 8, 3, 12, 14, 4, 9, 1, 6, 13],
                    [7, 1, 11, 8, 3, 6, 12, 9, 4, 13, 14, 5, 10, 2, 0],
                    [3, 8, 2, 12, 9, 4, 7, 13, 10, 5, 0, 14, 6, 11, 1],
                    [12, 4, 9, 3, 13, 10, 5, 8, 0, 11, 6, 1, 14, 7, 2],
                    [8, 13, 5, 10, 4, 0, 11, 6, 9, 1, 12, 7, 2, 14, 3],
                    [14, 9, 0, 6, 11, 5, 1, 12, 7, 10, 2, 13, 8, 3, 4],
                    [4, 14, 10, 1, 7, 12, 6, 2, 13, 8, 11, 3, 0, 9, 5],
                    [10, 5, 14, 11, 2, 8, 13, 7, 3, 0, 9, 12, 4, 1, 6],
                    [2, 11, 6, 14, 12, 3, 9, 0, 8, 4, 1, 10, 13, 5, 7],
                    [6, 3, 12, 7, 14, 13, 4, 10, 1, 9, 5, 2, 11, 0, 8],
                    [1, 7, 4, 13, 8, 14, 0, 5, 11, 2, 10, 6, 3, 12, 9],
                    [13, 2, 8, 5, 0, 9, 14, 1, 6, 12, 3, 11, 7, 4, 10],
                    [5, 0, 3, 9, 6, 1, 10, 14, 2, 7, 13, 4, 12, 8, 11],
                    [9, 6, 1, 4, 10, 7, 2, 11, 14, 3, 8, 0, 5, 13, 12],
                    [11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14]]
            print("Pbig", isLatin2(Pbig))
            print("Qbig", isLatin2(Qbig))
            print("Pbig,Qbig ortho", checkOrthogonality(Pbig, Qbig))
        Psmall, Qsmall = getPermutator(small)
        p = oddPlusOdd(Pbig, Psmall, False)
        if qrTest == 3 and engTest == "self":
            q = transposeMatrix(p)
            for i in range(3):
                for j in range(3):
                    q[big+i][big+j] = big + Qsmall[i][j]
        else:
            q = oddPlusOdd(Qbig, Qsmall, True)
        return p, q
        
    # Take care of hard-coded values of W that I just couldn't figure out any
    # decent way of handling otherwise.
    if W in hardcodedMOLS:
        P = hardcodedMOLS[W]["P"]
        if "Q" in hardcodedMOLS[W]:
            Q = hardcodedMOLS[W]["Q"]
        else:
            Q = transposeMatrix(P)
        return P, Q
        
    # Take care of some values of W for which I had to find good rules
    # essentially by trial and error.
    if W in [60, 84, 144, 156, 180, 192]:
        return getForDivisor(12)
    if W in [108]:
        return getForDivisor(27)
    if W in [132]:
        return getForDivisor(11)
    if W in [12, 36, 48, 162]:
        return getForDivisor(3)
    if W in [6]:
        return [], []
        
    # Now that the specific values of W are out of the way, fall back on some
    # general rules that seem to work pretty well.
    p = []
    q = []
    if W%2 == 1 and W%3 != 0:                 # W odd and not multiple of 3.
        p, q = oddNotDivisibleBy3(W)
    if p == [] and W%3 == 0 and W%6 != 0:     # 3|W but 6 doesn't divide W.   
        p, q = getForDivisor(3)
    if p == [] and pow2%2 == 0:               # W multiple of even power of 2.
        p, q = getForDivisor(4)
    if p == [] and pow2%2 == 1 and w == 1:    # W = odd power of 2.
        p, q = molsPower2(W)
    if p == [] and W % 18 == 0:
        p, q = getForDivisor(18)
        
    small = (W - 1) // 3
    big = W - small
    if p == [] and big == 2 * small + 1 and big % 2 !=0 and big % 3 != 0 \
            and small %2 != 0:
        Pbig, Qbig = getPermutator(big)
        Psmall, Qsmall = getPermutator(small)
        p = oddPlusOdd(Pbig, Psmall, False)
        q = oddPlusOdd(Qbig, Qsmall, True)
        
    # All else has failed.  Some last-ditch attempts.
    for divisor in [47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3]:
        if p == [] and W % divisor == 0:
            p, q = getForDivisor(divisor)
    if p != [] and q == []:
        q = transposeMatrix(p)
        
    return p,q

# Check if a square is latin.
def isLatin2(square):
    W = len(square)
    pure = list(range(W))
    for i1 in range(W):
        vector = []
        for i2 in range(W):
            vector.append(square[i1][i2])
        if pure != list(sorted(vector)):
            print("Row", i1, "is wrong:", vector)
            return False
    for i2 in range(W):
        vector = []
        for i1 in range(W):
            vector.append(square[i1][i2])
        if pure != list(sorted(vector)):
            print("Column", i2, "is wrong:", vector)
            return False
    return True

# Check if mutually orthogonal.
pairs = {}
def checkOrthogonality(square1, square2):
    global pairs
    okay = True
    pairs = {}
    W = len(square1)
    for i in range(W):
        for j in range(W):
            pair = str((square1[i][j], square2[i][j]))
            if pair in pairs:
                okay = False
                pairs[pair].append((i,j))
            else:
                pairs[pair] = [(i,j)]
    return okay

# Fetch the (pointer to) the lowest-level WxW plane of an N-cube, given the
# first N-2 coordinates.  Returns cube,i,plane, where cube is the WxWxW cube
# containing the target plane and i is the index within it of the plane.  If 
# index (which is a list) is empty, i.e. if N=2, then there is no containing 
# hypercube, and cube,i will simply be [],0.
#
# Notice, by the way, that if len(index)=N rather than N-2, then the returned
# "plane" will actually be the value of the cell at those coordinates rather
# than an actual WxW plane. 
def fetchPlane(hypercube, index):
    cube = []
    i = 0
    plane = hypercube
    for i in index:
        cube = plane
        plane = plane[i]
    return cube,i,plane

# Increment an index (which is a list of integers), used for a loop variable
# over the elements of a multi-dimensional array.  Returns either
# True (for end of loop) or False (not end of loop).  The index represents a
# list of coordinates, so "incrementing" it means incrementing the last one
# on the list and only proceeding to the earlier ones if the last one rolls 
# over.
def incrementIndex(index):
    if len(index) == 0:
        return False
    notDone = True
    for n in range(len(index)-1, -1, -1):
        index[n] += 1
        if index[n] < W:
            break
        index[n] = 0
        if n == 0:
            notDone = False
    return notDone

# Compare two multi-dimensional indexes -- i.e. arrays of numbers -- returning
# -1, 0, or +1, depending on whether the first index is "smaller than", 
# equal to, or "greater than" the second index.
def compareIndexes(index1, index2):
    for w in range(len(index1)):
        if index1[w] < index2[w]:
            return -1
        elif index1[w] > index2[w]:
            return +1
    return 0

# Check if an N-dimensional object m is magic or not.
def checkMagic(N, magic):
    error = False
    SN = (W * (W**N - 1)) / 2
    centerPoint = (W-1)/2
    # First check that every value in the range 0 to W**N-1 appears precisely
    # once.
    checkExcessiveMemory(W, N)
    values = bytearray(W**N) # 1 or 0 if value found or not found.
    index = [0]*N
    notDone = True
    count = 0
    while notDone:
        element = magic
        for m in index:
            element = element[m]
        if values[element] != 0:
            count += 1
            error = True
            print("Value", element, "duplicated", index)
            if count >= 10:
                print("...")
                break
        else:
            values[element] = 1
        # Update the loop variable, index[].
        notDone = incrementIndex(index)
    del values
    # Now check all stacks.  What this means is that for each
    # possible combination of N-1 coordinates held fixed, we let the Nth
    # coordinate vary. Once again this involves using index[] as a loop
    # variable.
    count = 0
    for fixedCoordinate in range(N):
        index = [0]*N
        notDone = True
        while notDone:
            if index[fixedCoordinate] == 0: 
                s = 0
                for i in range(W):
                    index[fixedCoordinate] = i
                    value = magic
                    for m in index:
                        value = value[m]
                    s += value
                index[fixedCoordinate] = "X"
                count += 1
                if s != SN:
                    print("Stack mismatch:", s, str(index).replace("'", ""))
                    error = True
                index[fixedCoordinate] = 0
            # Update the loop variable, index[].
            notDone = incrementIndex(index)
    print("%d stacks checked" % count)
    preDiagonalError = error
    if not onlyPartial:
        # Now check the principal diagonals.
        #
        # To understand how to determine the principal diagonals, consider the
        # N=3,W=3 case.  All principal diagonals must pass through the center 
        # of the hypercube, which is (1,1,1).  In addition to passing through
        # the center, it must also pass through 2 of the cells immediately
        # adjacent to the center cell, but on opposite sides of the center 
        # cell. For any given adjacent cell (x,y,z), the matching opposite 
        # cell is just (2-x,2-y,2-z). We can enumerate such adajacent cells 
        # easily (and in an automated way), since they're simply *all* of the 
        # cells in the 3x3x3 cube other than (1,1,1) itself.  So here are 
        # principal diagonals:
        #   (0,0,0) (1,1,1) (2,2,2)
        #   (0,0,1) (1,1,1) (2,2,1)
        #   (0,0,2) (1,1,1) (2,2,0)
        #   (0,1,0) (1,1,1) (2,1,2)
        #   (0,1,1) (1,1,1) (2,1,1)
        #   (0,1,2) (1,1,1) (2,1,0)
        #   (0,2,0) (1,1,1) (2,0,2)
        #   (0,2,1) (1,1,1) (2,0,1)
        #   (0,2,2) (1,1,1) (2,0,0)
        #   (1,0,0) (1,1,1) (1,2,2)
        #   (1,0,1) (1,1,1) (1,2,1)
        #   (1,0,2) (1,1,1) (1,2,0)
        #   (1,1,0) (1,1,1) (1,1,2)
        # Of course, 3 of these overlap with the rectilinear triplets already
        # checked.  I.e., they're not really diagonals.  But so what?  It 
        # doesn't hurt them to be checked again!
        #
        # How to extend this to N>3 (while retaining W=3) is obvious.  There 
        # will always be (3**N - 1)/2 diagonals, with N overlapping the 
        # rectilinears. Extending this to the case of odd W>2 (W=5, W=7, etc.) 
        # is quite easy as well, since the diagonals must be exactly those 
        # just listed, but extended by projecting and adding additional cells 
        # on each end.  In other words, there continue to be (3**N - 1)/2 - N 
        # diagonals, regardless of how big W may be.
        #
        # The case of even W is somewhat similar, but differs in detail. 
        # Consider the case of W=4.  For N=3, for example, the center of the 
        # hypercube is at (3/2, 3/2, 3/2), and the diagonals "pass through" 
        # this point though there's no actual cell there.  So to get the 
        # actual cells adjacent to the center point, we have to add all 
        # combinations of (±½, ±½, ±½) to the center coordinates. Actually, 
        # we'd be counding diagonals twice, since (for example) adding 
        # (-½, -½, -½) produces the same diagonal as adding
        # (+½, +½, +½) does.  Thus for N=3 there would be just 8/2=4 distinct 
        # diagonals, as opposed to (3**3-1)/2-3 = 10 for the W odd.
        #
        # And once again, the extension to N>3 for even W is straightforward:
        # 2**(N-1) distinct diagonals.
        
        # Step 1:  Find each of the principal diagonals, represented as an 
        # N-vector from the center point to an adjacent cell.  Store all of 
        # these N-vectors in a list called vectors[].
        zero = [0]*N
        odd = (W % 2 == 1)
        if odd:
            center = [(W-1)//2]*N
        else:
            center = [centerPoint]*N
        vectors = []
        if odd:
            # W is odd.
            # Construct vectors with individual ordinates having values -1, 0,
            # or +1, starting from [-1,...,-1], until reaching [0, ..., 0].
            vector = [-1]*N
            while True:
                vectors.append(copy.deepcopy(vector))
                for n in range(N-1, -1, -1):
                    vector[n] += 1
                    if vector[n] <= 1:
                        break
                    vector[n] = -1
                if compareIndexes(vector, zero) >= 0:
                    break
        else:
            # W is even.
            # Construct vectors with individual ordinates having values -1/2 
            # or +1/2, starting from [-1/2,...,-1/2], until surpassing 
            # [0, ..., 0].
            vector = [-1/2]*N
            while True:
                vectors.append(copy.deepcopy(vector))
                for n in range(N-1, -1, -1):
                    vector[n] += 1
                    if vector[n] <= 0.500000001:
                        break
                    vector[n] = -1/2
                if compareIndexes(vector, zero) >= 0:
                    break
        
        # Step 2:  Use the list of vectors, summing up the cell
        # values for each diagonal.
        for vector in vectors:
            # Use the "vector" to create a list of coordinates (diagonal[]) 
            # for the cells comprising the diagonal
            if odd: # W odd
                diagonal = [center]
                for w in range((W-1)//2):
                    cell = []
                    for n in range(N):
                        cell.append(center[n] + (w+1)*vector[n])
                    diagonal.append(cell)
                    cell = []
                    for n in range(N):
                        cell.append(center[n] - (w+1)*vector[n])
                    diagonal.append(cell)
            else: # W even
                diagonal = []
                for w in range(W//2):
                    cell = []
                    for n in range(N):
                        cell.append(round(center[n] + (2*w+1)*vector[n]))
                    diagonal.append(cell)
                    cell = []
                    for n in range(N):
                        cell.append(round(center[n] - (2*w+1)*vector[n]))
                    diagonal.append(cell)
            s = 0
            for index in diagonal:
                cube,i,value = fetchPlane(magic, index)
                s += value
            if s != SN:
                def key(element):
                    sum = 0
                    for e in element:
                        sum = sum*W + e
                    return sum
                diagonal.sort(key=key)
                print("Diagonal mismatch:", s, str(diagonal).replace("'", ""))
                error = True
        print("%d diagonals checked" % len(vectors))
        del vectors
    return error, preDiagonalError
    
if testMOLS:
    print("Testing MOLS q+r construction only.")
    for W in range(3, 1001):
        p, q = getPermutator(W)
        if p != []:
            print("W =", W, "already implemented.")
            continue
        for big in smallPrimes:
            if big < 3:
                continue
            small = W - big
            if small < 3:
                break
            if small > big:
                continue
            msg = "W = %d + %d = %d: " % (big, small, W)
            Pbig, Qbig = getPermutator(big)
            if Pbig == []:
                msg += " Big block unimplemented."
                print(msg)
                continue
            Psmall, Qsmall = getPermutator(small)
            if Psmall == []:
                msg += " Small block unimplemented."
                print(msg)
                continue
            p = oddPlusOdd(Pbig, Psmall, False)
            q = oddPlusOdd(Qbig, Qsmall, True)
            if isLatin2(p):
                msg += " P is latin."
            else:
                msg += " P is not latin."
            if isLatin2(q):
                msg += " Q is latin."
            else:
                msg += " Q is not latin."
            if checkOrthogonality(p, q):
                msg += " Mutually orthogonal."
            else:
                msg += " Not mutually orthogonal."
            print(msg)   
    sys.exit(0)

p,q = getPermutator(W)
if p == []:
    printReport("Not yet implemented", 2)
    sys.exit(1)

if debug:
    def isSelfComplementary(p):
        W = len(p)
        for i in range(W):
            for j in range(W):
                if p[i][j] != W - 1 - p[W - 1 - i][W - 1 - j]:
                    return False
        return True
    
    if not isLatin2(p):
        print("p is not latin")
    if not isLatin2(q):
        print("q is not latin")
    if not checkOrthogonality(p,q):
        print("p and q not MOLS")
        if False:
            print ("Duplicates")
            for i in sorted(pairs):
                if len(pairs[i]) > 1:
                    print("\t%s\t" % i, pairs[i])
        else:
            # For a prettier picture.
            print ("Failures (. = okay, X = duplicate):")
            pairsSquare = []
            byRow = []
            missing = []
            for i in range(W):
                pairsSquare.append(['.']*W)
                byRow.append([])
            for i in range(W):
                for j in range(W):
                    key = str((i,j))
                    if key not in pairs:
                        missing.append(key)
            for pair in sorted(pairs):
                if len(pairs[pair]) > 1:
                    pairSymbol = "X"
                    for position in pairs[pair]:
                        pairsSquare[position[0]][position[1]] = pairSymbol
                        byRow[position[0]].append((position[1], pair))
            for i in range(4):
                msg = "      "
                for j in range(W):
                    msg += " " + ("%4d" % j)[i]
                print(msg)
            msg = "      "
            for j in range(W):
                msg += " ="
            print(msg)
            for i in range(W):
                row = pairsSquare[i]
                msg = "%4d: " % i
                for element in row:
                    msg += " " + element
                print(msg)
            for i in range(W):
                print("Row =", i, 
                      ", Duplicates (col,P,Q) =", sorted(byRow[i]))
            print("Missing pairs:", missing)
    pSumDiagonal = 0
    pSumAntidiagonal = 0
    qSumDiagonal = 0
    qSumAntidiagonal = 0
    pDiagonal = []
    pAntidiagonal = []
    qDiagonal = []
    qAntidiagonal = []
    for i in range(W):
        pDiagonal.append(p[i][i])
        pSumDiagonal += p[i][i]
        pAntidiagonal.append(p[i][W-1-i])
        pSumAntidiagonal += p[i][W-1-i]
        qDiagonal.append(q[i][i])
        qSumDiagonal += q[i][i]
        qAntidiagonal.append(q[i][W-1-i])
        qSumAntidiagonal += q[i][W-1-i]
    pDiagonal.sort()
    pAntidiagonal.sort()
    qDiagonal.sort()
    qAntidiagonal.sort()
    targetSum = (W*(W-1))//2
    print("p,q target sum", targetSum)
    print("p diagonal", pSumDiagonal, pDiagonal)
    print("p antidiagonal", pSumAntidiagonal, pAntidiagonal)
    print("q diagonal", qSumDiagonal, qDiagonal)
    print("q antidiagonal", qSumAntidiagonal, qAntidiagonal)
    print(". Width", W)
    if sorted(pDiagonal) == list(range(W)) \
            and sorted(pAntidiagonal) == list(range(W)):
        print(". P diagonal LS")
    else:
        print(". P not diagonal")
    if sorted(qDiagonal) == list(range(W)) \
            and sorted(qAntidiagonal) == list(range(W)):
        print(". Q diagonal LS")
    else:
        print(". Q not diagonal")
    if pSumDiagonal == targetSum and pSumAntidiagonal == targetSum:
        print(". P sums okay")
    else:
        print(". P sums bad")
    if qSumDiagonal == targetSum and qSumAntidiagonal == targetSum:
        print(". Q sums okay")
    else:
        print(". Q sums bad")
    if p == transposeMatrix(q):
        print(". Self-orthogonal")
    else:
        print(". Not self-orthogonal")
    if isSelfComplementary(p):
        print(". P self-complementary")
    else:
        print(". P not self-complementary")
    if isSelfComplementary(q):
        print(". Q self-complementary")
    else:
        print(". Q not self-complementary")
    print()
    print("p = [", end="")
    for i in range(len(p)):
        end = ","
        if i == len(p)-1: end = "]"
        print("\t", p[i], end)
    print("q = [", end="")
    for i in range(len(q)):
        end = ","
        if i == len(q)-1: end = "]"
        print("\t", q[i], end)

# M is a semi-magic hypercube of some dimension N and width W.
# A is a latin hypercube of the same dimension N and width W as M.
# P and Q are MOLS of width W. Computes semi-magic M' and latin A' 
# of the same width W but incremented dimension N+1.
def constructDimension(M, A, P, Q):

    # Recursively combine 2 multi-dimensional arrays E and A
    # identical dimension and width as E = E + WtoN * A.
    # E is modified in place, but also returned.
    def combine(E, WtoN, A):
        for i in range(len(E)):
            if type(E[i]) == list:
                combine(E[i], WtoN, A[i])
            else:
                E[i] += WtoN * A[i]
        return E

    # Find W and compute W**N.
    W = len(P)
    WtoN = 1
    v = M
    while type(v) == list:
        WtoN *= len(v)
        v = v[0]
        
    # Construct A' from A.
    Aprime = []
    for i in range(W):
        row = []
        for j in range(W):
            row.append(A[Q[i][j]])
        Aprime.append(row)
        
    # Construct M' from M and A'.
    E = []
    for i in range(W):
        row = []
        for j in range(W):
            row.append(copy.deepcopy(M[P[i][j]]))
        E.append(row)
    Mprime = combine(E, WtoN, Aprime)
    
    return Mprime, Aprime

# Main loop =============================================================
# Construct higher-dimensional objects for dimensions 2 through maxN, 
# starting from the given seeds for dimension 1.
wPower = 1
for N in range(2, maxN+1):
    startTime = time.time()
    print()
    print("Dimension", N)
    print("Target sum", (W*(W**N-1))//2)
    predictExcessiveMemory(W, N)
    m, a = constructDimension(m, a, p, q)
    # Check the result.
    error, rectError = checkMagic(N, m)
    if error and rectError:
        printReport("Failure", N)
    elif onlyPartial and not error:
        print(".", N, "Semi-magic")
        printReport("Semi-magic success", N)
    elif error and not rectError:
        print(".", N, "Semi-magic")
        printReport("Semi-magic success, magic failure", N)
    else:
        print(".", N, "Magic")
        printReport("Magic success", N)

    # Print the results
    if not quiet:
        if python:
            print("m =", m)
        else:
            # Print a multi-dimensional object.
            def printMultiDimensional(m):
                for s in m:
                    if type(s[0]) != list:
                        for e in s:
                            print("%5d" % e, end="")
                        print()
                    else:
                        printMultiDimensional(s)
                        print()
        
            print("M")
            printMultiDimensional(m)  
            if debug:
                print()
                print("A")
                printMultiDimensional(a)
