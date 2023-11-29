from math import sqrt, inf
from sympy import *
import numpy as np

global x
x = Symbol("x")

def calculateAll(coeffs, values):
    out = []
    for v in values:
        out += [(v,calculate(coeffs,v))]
    return out

def calculate(coeffs, x):
    out = 0
    for index, i in enumerate(coeffs):
        if index > 0:
            num = i * (x ** index)
            out += num
    out += coeffs[0]
    return out

def mkTerm(coeff, deg): 
    if coeff == 0 or deg < 0:
        return
    coeff = abs(coeff)
    x = ""
    if deg >= 1:
        x = "x"
    out = ""
    if (x == "x" and coeff != 1) or x == "":
        out += str(round(coeff, 2))
        out += "*"
    out += x
    if deg > 1:
        out += "**" + str(deg)
    out = list(out)
    if out[len(out)-1] == "*":
        out[len(out)-1] = ""
    return "".join(out)

def formatPoly(coeffs):
    out = ""
    outList = []

    for i in range(len(coeffs)): # for each coefficient formatPoly is passed
        if coeffs[i] < 0:
            t = mkTerm(coeffs[i], i)
            if t is not None:
                outList.append("-" + t)
        else:
            t = mkTerm(coeffs[i], i)
            if t is not None:
                outList.append(t)

    outList.reverse() # we added everything in by ascending order of degree, now reverse it
   
    for i in range(len(outList)):
        if outList[i] is not None:
            out += str(outList[i])
            if i < len(outList) - 1:
                out += " + "

    out = list(out)
    for i in range(len(out)):
        if out[i] == "+" and out[i+1] == " " and out[i+2] == "-": # if there's somewhere that a plus sign is in front of a number's negative sign
            # alter the terms so that only the negative sign exists            
            out[i] = "-"
            out[i+2] = ""
            
    return "".join(out)

def getDerivative(coeffs):
    f = formatPoly(coeffs)
    return diff(f, x)
    
def newtonsMethod(tolerance, function, initialGuess):
    '''
    An intelligent method of performing a guess-and-check. Follows the formula:

    x_(n+1) = x_n - (f(x_n)/f'(x_n))

    where f(x) is the given function, f'(x) is its derivative, n is any step number and n+1 is the next step.
    Uses the SymPy library to compute derivatives and create functions out of our expressions.
    '''
    guess = initialGuess # uses some arbitrary number to create a starting guess
    prevGuess = initialGuess
    f_xn = lambdify(x, formatPoly(coeffs)) # create a function that takes an x-value and returns f(x), computed from "coeffs"
    fprime_xn = lambdify(x, getDerivative(coeffs)) # same thing but it returns f'(x)
    i = 0
    maxIterations = 500 # if it diverges away from our zero, we don't want it to go on for forever
    while i < maxIterations:
        numerator = f_xn(guess) # creates f(x_n)
        denominator = fprime_xn(guess) # creates f'(x_n)
        if denominator != 0: # if we're not dividing by zero
            guess = guess - (numerator / denominator)
        else: # if we are trying to divide by zero, something must have gone wrong
            return None
        if abs(calculate(coeffs, guess)) < tolerance: # if we're within our tolerance of accuracy, consider the answer to be correct
            return guess # and return it
        prevGuess = guess
        i += 1
    return None
    
def zeros(coeffs):
    out = []
    if len(coeffs) ==  1: # degree 0 polynomial
        if coeffs[0] == 0:
            out.append(inf)

    elif len(coeffs) == 2: # degree 1 polynomial; use general equation
        out.append(-coeffs[0] / coeffs[1])

    elif len(coeffs) == 3: # degree 2 polynomial; use quadratic solution
        a = coeffs[2]
        b = coeffs[1]
        c = coeffs[0]
        minusB = -b
        num = (b**2) - (4 * a * c)
        if num > 0:
            num = sqrt(num)
            minusB /= 2 * a
            num /= 2 * a
            out.append(minusB - num)
            out.append(minusB + num)
        elif num == 0:
            out.append((-b / (2 * a)))
    elif len(coeffs) > 3: # works for any n-degree polynomial; use Newton's method of approximation
        maxRange = 20
        guessPrecision = 0.5
        guesses = np.arange(-abs(maxRange), abs(maxRange) + 1, guessPrecision) # controls the range of our guesses and how often our guesses are
        for i in guesses:
            zero = newtonsMethod(tolerance, formatPoly(coeffs), i)
            if zero is not None:
                zero = round(zero, 8)
                if zero == 0:
                    zero = abs(zero) # ensure that we don't give "-0.0" as an answer just because it looks better
                if zero not in out: # if we don't already have it
                    out.append(zero) # add it to our list
        out.sort()
    else:
        pass
    return out


tolerance = 1e-10
coeffs = [8, -20, 14, -6, 1]

print(f"{zeros(coeffs)} is/are the location of the zeros.")
