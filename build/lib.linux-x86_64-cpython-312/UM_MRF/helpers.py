from numpy import abs

def isapprox(a, b, atol=1e-9):
    return (abs(a - b) <= atol)

def isnapprox(a, b, atol=1e-9):
    return (abs(a - b) > atol)