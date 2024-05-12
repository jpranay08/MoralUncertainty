import itertools
import time
import random
def safe_divide(numerator, denominator):
    if (denominator < 0) | (denominator == 0):
        return numerator
    else:
        return (numerator / denominator)


# default prior values are uniform uninformative
def massComb(masses, prior0=0.5, prior1=0.5, prior01=0):

    # the space in the hypothesis space that is not in the evidence space
    prior_theta = 1 - (prior0 + prior1 + prior01)

    # since sum of m(A) must equal 1 there may be frame of descernment is
    # what's left over
    for i in range(len(masses)):
        masses[i].append( 1 - sum(masses[i]))
    

    ##########################
    # PERFORM MASS COMBINATION
    #########################
    #print("The different mass functions are:")
    intrsxn_array_dim = len(masses[1])
    # set the dimenensions of the mass comb. matrix.
    intrsxn_array = [
        [0 for j in range(intrsxn_array_dim)] for i in range(intrsxn_array_dim)]

    ############### BEGIN: Combining all bpa's ##############################
    for i in range(0,len(masses)-1):
        if i == 0:
            K = 1  
            m0 = masses[0][0]
            m1 = masses[0][1]
            m01 = masses[0][2]
            m_theta = masses[0][3]
            #print("First mass assignment. m0:", m0, "m1: ", m1,
            #      "m01: ", m01, "m_theta: ", m_theta, "K: ", K)

        new_mass = [[m0, m1, m01, m_theta]]
    
        for col in range(intrsxn_array_dim):
            for row in range(intrsxn_array_dim):
                intrsxn_array[row][col] = round(new_mass[0][col]*masses[i + 1][row],2)
        #print(new_mass,"new mass")
        

        # CALCULATE K - the measure of conflict
        K = intrsxn_array[0][1] + intrsxn_array[1][0]
        if(K==1.0):
            K=0.9
        # Calculate belief functions
        m0 = (intrsxn_array[0][0] + intrsxn_array[2][0] + intrsxn_array[3][0]
            + intrsxn_array[0][2] + intrsxn_array[0][3]) / (1 - K)
        m1 = (intrsxn_array[1][1] + intrsxn_array[1][2] + intrsxn_array[1][3]
            + intrsxn_array[2][1] + intrsxn_array[3][1]) / (1 - K)
        m01 = intrsxn_array[2][3] / (1 - K)
        # normalize to emphasise agreement
        m_theta = intrsxn_array[3][3] / (1 - K)
        #print("Next mass assignment. m0:", m0, "m1: ", m1,
        #      "m01: ", m01, "m_theta: ", m_theta, "K: ", K)
    ############### END: Combining all bpa's ###############################

    #print("m0:", m0, "m1: ", m1, "m01: ",
    #      m01, "m_theta: ", m_theta, "K: ", K)
    # INCLUDE PRIOR INFORMATION
    #print("\n")
    #print("prior0: ", prior0, "prior1: ", prior1,
    #      "prior01: ", prior01, "prior_theta: ", prior_theta)
    # basic certainty assignment (bca) and normalize
    certainty_denominator = (safe_divide(numerator = m0, denominator = prior0) + safe_divide(numerator = m1, denominator = prior1)
                             + 
                             safe_divide(
                                 numerator=m01, denominator=prior01)
                             + safe_divide(numerator=m_theta, denominator=prior_theta))
    # print(certainty_denominator)
    # time.sleep(2)
    if(certainty_denominator==0):
        certainty_denominator=0.5
    C0 = safe_divide(numerator=m0, denominator=prior0) / \
        certainty_denominator
    C1 = safe_divide(numerator=m1, denominator=prior1) / \
        certainty_denominator
    C01 = safe_divide(
        numerator=m01, denominator=prior01) / certainty_denominator
    C_theta = safe_divide(
        numerator=m_theta, denominator=prior_theta) / certainty_denominator
    #print("C0: ", C0, "C1: ", C1, "C01: ", C01, "C_theta: ", C_theta)
    #print("C0 + C1 + C01 + C_theta: ", C0 + C1 + C01 + C_theta, "\n")

    #print("inrsxn_array:", intrsxn_array[0][1])
    
    blf0 = round(m0,2)
    blf1 = round(m1,2)
    blf01 = round(m0 + m1 + m01,2)

    plsb0 = round(m0 + m01 + m_theta,2)
    plsb1 = round(m1 + m01 + m_theta,2)
    plsb_theta = 1

    mass_fxn_values = {"blf0": blf0, "blf1": blf1, "blf01": blf01, \
                        "plsb0": plsb0, "plsb1": plsb1, "plsb_theta": plsb_theta}

    return mass_fxn_values


def tMassComb(masses, prior_a=0.33, prior_b=0.33, prior_c=0.33, prior_ab=0, prior_bc=0, prior_ac=0, prior_abc=0):
    # the space in the hypothesis space that is not in the evidence space
    prior_abc = 1 - (prior_a + prior_b + prior_c + prior_ab + prior_bc + prior_ac)

    # since sum of m(A) must equal 1 there may be frame of discernment is
    # what's left over
    for i in range(len(masses)):
        masses[i].append(1 - sum(masses[i]))

    ##########################
    # PERFORM MASS COMBINATION
    #########################

    intrsxn_array_dim = len(masses[1])
    # set the dimensions of the mass comb. matrix.
    intrsxn_array = [
        [0 for j in range(intrsxn_array_dim)] for i in range(intrsxn_array_dim)]
    for i in range(0, len(masses) - 1):
        if i == 0:
            K = 1
            m_a = masses[0][0]
            m_b = masses[0][1]
            m_c = masses[0][2]
            m_ab = masses[0][3]
            m_bc = masses[0][4]
            m_ac = masses[0][5]
            m_abc = masses[0][6]

        new_mass = [[m_a, m_b, m_c, m_ab, m_bc, m_ac, m_abc]]
        for col in range(intrsxn_array_dim):
            for row in range(intrsxn_array_dim):
                intrsxn_array[row][col] = round(
                    new_mass[0][col] * masses[i + 1][row], 2)

        # CALCULATE K - the measure of conflict
        K = 0
        for col in range(intrsxn_array_dim):
            for row in range(intrsxn_array_dim):
                if col != row:
                    K += intrsxn_array[row][col]

        # Calculate belief functions
        m_a = sum(intrsxn_array[0])
        m_b = sum(intrsxn_array[1])
        m_c = sum(intrsxn_array[2])
        m_ab = intrsxn_array[3][0] + intrsxn_array[3][1] + intrsxn_array[3][2]
        m_bc = intrsxn_array[4][0] + intrsxn_array[4][1] + intrsxn_array[4][2]
        m_ac = intrsxn_array[5][0] + intrsxn_array[5][1] + intrsxn_array[5][2]
        m_abc = intrsxn_array[6][0] + intrsxn_array[6][1] + intrsxn_array[6][2]

    # INCLUDE PRIOR INFORMATION
    # basic certainty assignment (bca) and normalize
    certainty_denominator = (
        safe_divide(m_a, prior_a)
        + safe_divide(m_b, prior_b)
        + safe_divide(m_c, prior_c)
        + safe_divide(m_ab, prior_ab)
        + safe_divide(m_bc, prior_bc)
        + safe_divide(m_ac, prior_ac)
        + safe_divide(m_abc, prior_abc)
    )

    if certainty_denominator == 0:
        certainty_denominator = 0.5

    C_a = safe_divide(m_a, prior_a) / certainty_denominator
    C_b = safe_divide(m_b, prior_b) / certainty_denominator
    C_c = safe_divide(m_c, prior_c) / certainty_denominator
    C_ab = safe_divide(m_ab, prior_ab) / certainty_denominator
    C_bc = safe_divide(m_bc, prior_bc) / certainty_denominator
    C_ac = safe_divide(m_ac, prior_ac) / certainty_denominator
    C_abc = safe_divide(m_abc, prior_abc) / certainty_denominator

    mass_fxn_values =[round(C_a,2),round(C_b,2),round(C_c,2)]
    return mass_fxn_values


def makeMf(a,i):
  res=[0,0,0,0,0,0]
  res[i]=a
  b=round(random.uniform(0,1-a),2)
  c=round(1-a-b,2)
  if(i==0):
    res[1]=b
    res[2]=c
  elif(i==1):
    res[0]=b
    res[2]=c
  elif(i==2):
    res[0]=b
    res[1]=c
  return res


def getDempsterCred(arr):
  imass=[]
  for i in range(len(arr)):
    ma=makeMf(arr[i],i)
    imass.append(ma)

  return tMassComb(imass)



