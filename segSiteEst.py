__author__ = 'ruhuajiang'
'''
Using simple monte carlo to estimate segregating sites.
'''
import matplotlib.pyplot as plt
import numpy as np
import math


def get_prob(theta,k,t_total):
    return (math.pow((float(theta*t_total)/2),k))*(math.exp(-theta*t_total/2))/(math.factorial(k))


def monte_carlo(theta,k,sample_size, exp_count):
    #Repeat exp_count number of times
    p_list = []
    for _ in range(exp_count):
        t_total = 0
        for i in range(sample_size,1,-1):
            scale = float(2)/(i*(i-1))
            t = np.random.exponential(scale,1)
            t_total += i*t
        #print "t_total",t_total
        p = get_prob(theta, k, t_total)
        p_list.append(p)
    return float(sum(p_list))/len(p_list)
    #Get the histogram
    #plt.hist(p_list)
    #plt.show()

def main():
    print "Genealogy simulation starts ..."
    sample_size = 10
    theta = 1
    k = 2
    x_arr = range(10,10000,100)
    y_arr = []
    y2 = []
    for exp_count in x_arr:
        prob = monte_carlo(theta,k,sample_size,exp_count)
        y_arr.append(prob)
        y2.append(0.2135)
    plt.xlabel("R, number of experiments")
    plt.ylabel("Probability, P{S = 2; theta = 1}")
    plt.plot(x_arr,y_arr)
    plt.plot(x_arr,y2,'r')
    plt.show()


if __name__ == "__main__":
    main()