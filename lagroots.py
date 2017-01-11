from numpy.polynomial.laguerre import lagroots
import csv

with open("C:\\Users\\User\\Documents\\Cpp\\integration\\lagroots.csv", mode='w') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 'root'])
    for i in range(1, 81):
        coef = [0]*(i+1)
        coef[i] = 1
        for r in lagroots(coef):
            writer.writerow([i,"{0:.40f}".format(r)])
