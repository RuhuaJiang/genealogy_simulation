__author__ = 'ruhuajiang'
'''
This program simulate genealogy based on the sequence data
'''
import numpy as np
import random

debug_fxn_list = []
def process_seq():
    '''Get the gene tree'''
    X = [
       [4,0],
       [3,4,0],
       [2,1,3,4,0],
       [0]
    ]
    n = [2, 1,1,1]
    return X,n

def isDistinct(hap,X,i):
    lst = [e for e in X if e != hap]
    other_mutations = set()
    for l in lst:
        for x in l:
            if x!=0:
                other_mutations.add(x)
    #print "hap",hap,"other mutations",other_mutations
    if hap[0] in other_mutations:
        return False
    else:
        return True

def becomeAnother(hap,X,n):
    try:
        idx = X.index(hap[1:])
        return idx
    except ValueError:
        return -1

def scan(X,n):
    coalescent_criteria = []
    no_count_change_mut= []
    k_j_pair_mut = []
    for i in range(len(X)):
        if X[i] == [0]:
            continue
        if n[i] >= 2:
            coalescent_criteria.append(i)
        elif n[i] == 1:
            if isDistinct(X[i],X,i):
               j = becomeAnother(X[i],X,n)
               if j != -1:
                   k_j_pair_mut.append((i,j))
               else:
                   no_count_change_mut.append(i)
        else:
            print "ERROR"
            exit()
    return coalescent_criteria,no_count_change_mut,k_j_pair_mut

def sample_prob_list(prob_list):
    #print "Transit probability list:", prob_list
    if len(prob_list) ==0:
       return -1
    idx = 0
    for i in range(1,len(prob_list)):
        prob_list[i] += prob_list[i-1]
    rd = random.random()
    if rd < prob_list[0]:
        idx = 0
    for i in range(1,len(prob_list)):
        if rd >= prob_list[i-1] and rd <= prob_list[i]:
            idx = i
    return idx

def single_chain(X,n,theta):
    #print "\n======"
    #print "X",X
    #print "n",n
    N = sum(n)
    #print "N",N
    if X == [[0]] and n == [2]:
        return 1

    coalescent_criteria,no_count_change_mut,k_j_pair_mut = scan(X,n)
    #print "Index of coalesce, mutate1 , mutate2: ",coalescent_criteria,no_count_change_mut,k_j_pair_mut
    first_term = 0
    for k in coalescent_criteria:
        first_term += float(n[k]*(n[k]-1))/(N*(N-1 + theta))
    second_term = len(no_count_change_mut)*float(theta)/(N*(N-1+theta))
    third_term = float(len(k_j_pair_mut))*theta/(N*(theta + N -1))
    fxn = first_term + second_term + third_term
    #print "fxn:",fxn
    #debug_fxn_list.append(fxn)
    #print "DEBUG_FXN_LIST:",debug_fxn_list
    #print "Three probability terms:",first_term,second_term,third_term
    prob_list = []
    for k in coalescent_criteria:
        prob_list.append(float(n[k]*(n[k]-1))/(fxn*N*(theta+N-1)))
    for i in range(len(no_count_change_mut)+len(k_j_pair_mut)):
        prob_list.append(float(theta)/(fxn*N*(theta+N-1)))
    count_list = [len(coalescent_criteria),len(no_count_change_mut),len(k_j_pair_mut)]
    idx = sample_prob_list(prob_list)
    for i in range(1,len(count_list)):
        count_list[i] += count_list[i-1]
    if idx < count_list[0]:
        #coalescent event
        #print "Event:coalesce"
        k = coalescent_criteria[idx]
        n[k] -= 1
        return fxn*single_chain(X,n,theta)

    elif idx >= count_list[0] and idx < count_list[1]:
        #first kind mutation event
        k = no_count_change_mut[idx - count_list[0]]
        #print "Event:first kind mutate"
        X[k] = X[k][1:]
        return fxn*single_chain(X,n,theta)
    else:
        #second kind mutation event
        k,j = k_j_pair_mut[idx - count_list[1]]
        #print "Event:second kind mutate"
        X.pop(k)
        n[j] +=1
        n.pop(k)
        return fxn*single_chain(X,n,theta)

def important_sampling(X,n,r):
    print "Method: important sampling."
    theta = 1
    prob = 0
    for  i in range(r):
      X_dup = X[:]
      n_dup = n[:]
      p =single_chain(X_dup,n_dup,theta)
      #print p
      prob += p
    prob = float(prob)/r
    print "Number of replicates:", r
    print "Probability is:",prob


def main():
    print "Genealogy simulation starts ..."
    X,n  = process_seq()
    r = 10000
    #monte_carlo(data)
    important_sampling(X,n, r)

if __name__ == "__main__":
    main()