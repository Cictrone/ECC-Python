class ECC():

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p
        self.r_x = None
        self.r_y = None

    def setPoint(self, z):
        self.g_x = z[0]
        self.g_y = z[1]

    def intModInverse(self, e):
        return pow(e, self.p-2, self.p)

    def PointDecompress(self, P):
        x = P[0]
        if (self.p % 4 != 3):
            print("Square root not implemented for prime 1 mod 4")
        z = x**3 + self.a*x + self.b
        z %= self.p
        y = pow(z, (self.p+1)//4, self.p)
        i = P[1]
        y_i = y % 2
        if (i==y_i):
            return (x,y)
        return (x, self.p-y)

    def raiseByExponent(self, exp):
        if exp is 1:
          self.r_x = self.g_x
          self.r_y = self.g_y
          return
        gArray = []
        gArray.append(self.g_x)
        gArray.append(self.g_y)
        self.add(self.g_x, self.g_y, self.g_x, self.g_y)
        gArray.append(self.r_x)
        gArray.append(self.r_y)

        for i in range(2, exp.bit_length()):
            self.add(self.r_x, self.r_y, self.r_x, self.r_y)
            gArray.append(self.r_x)
            gArray.append(self.r_y)

        lowestBit = 0
        for i in range(exp.bit_length()):
            if ((exp & (1 << i)) != 0):
              lowestBit = i
              break
        self.r_x = gArray[2*lowestBit]
        self.r_y = gArray[(2*lowestBit)+1]
        lowestBit += 1
        for i in range(lowestBit, exp.bit_length()):
            if ((exp & (1 << i)) != 0):
              self.add(self.r_x, self.r_y, gArray[2*i], gArray[(2*i)+1])

    def add(self, x_0, y_0, x_1, y_1):
        if x_0 is None:
          self.r_x = x_1
          self.r_y = y_1
          return
        if x_1 is None:
          self.r_x = x_0
          self.r_y = y_0
          return
        if (x_0 == x_1) and ((y_1 + y_0) == self.p):
          self.r_x = None
          self.r_y = None
          return
        if x_0 is not x_1:
          lambduh = ((y_1 - y_0) % self.p)
          lambduh *= self.intModInverse(x_1 - x_0)
          lambduh %= self.p
          self.r_x = (lambduh*lambduh - x_0 - x_1) % self.p
          self.r_y = (lambduh*(x_0 - self.r_x) - y_0) % self.p
          return
        lambduh = (3*x_1*x_1 + self.a) % self.p
        lambduh *= self.intModInverse(y_0*2) % self.p
        self.r_x = (lambduh*lambduh - x_0 - x_1) % self.p
        self.r_y = ((x_0 - self.r_x)*lambduh - y_0) % self.p


if __name__ == '__main__':
    curve = ECC(4, 20, 29)
    #P = curve.PointDecompress((18,1))
    #print(P)
    curve.setPoint((8,10))
    curve.raiseByExponent(9)
    print(curve.r_x, curve.r_y)
