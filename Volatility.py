import csv, math
import matplotlib.pyplot as plt
import numpy as np

class parseVol:
    x_values = []
    y_values = []
    contract = None
    expiry = None
    def __init__(self, filename):
        self.filename = filename

    def parsefile(self):
        with open(self.filename, newline='') as csvfile:
            reader =csv.reader(csvfile)
            rownumber = 0
            for row in reader:
                rownumber = rownumber + 1
                if rownumber == 1:
                    for i in range(2,len(row)):
                        self.x_values.append(row[i])
                elif rownumber == 2:
                    self.contract = row[0]
                    self.expiry = row[1]
                    for i in range( 2, len( row ) ):
                        self.y_values.append(row[i])

    def get_expiry(self):
        return self.expiry

    def handle_nan(self, value):
        if  value == '':
            return 0
        else:
            return value

    def initialize(self):
        self.xs = []
        self.ys = []
        n = len( self.x_values )

        self.yt = []
        u = []

        for i in range(n + 1):
            self.xs.append(0)
            self.ys.append(0)

        for i in range(n + 1):
            self.yt.append(0)

        for i in range(n):
            u.append(0)

        self.xs[0] = 0
        self.ys[0] = 0

        for j in range( 1, n + 1):
            self.xs[j] = float(self.x_values[j - 1])
            self.ys[j] = float(self.handle_nan(self.y_values[j - 1]))

        if self.xs[0] < self.xs[1]:
            while (abs( self.xs[n]) < 1.0e-15 and n > 1):
                n = n - 1

        p = qn = sig = un = 0.0
        u[0] = u[1] = self.yt[0] = self.yt[1] = 0.0

        for i in range( 2, n ):
            sig = (self.xs[i] - self.xs[i - 1]) / (self.xs[i + 1] - self.xs[i - 1])
            p = sig * self.yt[i - 1] + 2.0
            self.yt[i] = (sig - 1.0) / p
            d = (self.ys[i + 1] - self.ys[i]) / (self.xs[i + 1] - self.xs[i]) - (self.ys[i] - self.ys[i - 1]) / (self.xs[i] - self.xs[i - 1])
            u[i] = (6.0 * d / (self.xs[i + 1] - self.xs[i - 1]) - sig * u[i - 1]) / p

        self.yt[n] = (un - qn * u[n - 1]) / (qn * self.yt[n - 1] + 1.0)

        for i in range( n - 1, 0, -1 ):
            self.yt[i] = self.yt[i] * self.yt[i + 1] + u[i]

    def interpolate(self, x_value):
    #cubic spline#
        n = len( self.x_values )
        if self.xs[1] == self.xs[n]:
            raise Exception('xs[1]==xs[n]')
        elif self.xs[1] < self.xs[n]:
            if x_value <= self.xs[1]:
                return self.ys[1]
            elif x_value >= self.xs[n]:
                return self.ys[n]
        else:
            if x_value >= self.xs[1]:
                return self.ys[1]
            elif x <= self.xs[n]:
                return self.ys[n]

        klo = 1
        khi = n
        k = khi - klo
        while k > 1:
            k = khi - klo
            if self.xs[k] > x_value:
                khi = k
            else:
                klo = k
                k = khi - klo

        h = self.xs[khi] - self.xs[klo]
        a = (self.xs[khi] - x_value) / h
        b = (x_value - self.xs[klo]) / h
        d = a * self.ys[klo] + b * self.ys[khi] + ((math.pow( a, 3 ) - a) * self.yt[klo] + (math.pow( b, 3 ) - b) * self.yt[khi]) * math.pow( h,
                                                                                                                      2 ) / 6.0
        return d

    def interpolate_volatilities(self):
        xpoints = []
        xpoints.append(float(self.x_values[0]))
        maximum = float(self.x_values[len(self.x_values)-1])
        #print(maximum)
        x= xpoints[0]
        while (x + 0.1 <= maximum):
            x = x + 0.1
            xpoints.append(x)
        ypoints = []
        for x in xpoints:
            ypoints.append(self.interpolate(x))

        with open( 'output.csv', 'w' ) as csvfile:
            writer = csv.writer( csvfile)
            writer.writerow(xpoints)
            writer.writerow(ypoints)

    def printstuffs(self):
        print(self.x_values)
        print(self.y_values)
        #for x in range(3, 0, -1):
        #    print(x)

def main():
    print("this is the main")
    myobj = parseVol('noisy.csv')
    myobj.parsefile()
    myobj.initialize()
    myobj.printstuffs()
    print(myobj.interpolate(599.8))
    myobj.interpolate_volatilities()

#if __name__=="__main__":
 #   main()