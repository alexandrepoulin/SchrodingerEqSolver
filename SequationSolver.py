##imports for plotting
import mpmath as mp
import matplotlib.pyplot as plt
import numpy as np

##global variables for the program

STEPSIZE    = 0.0001        ##the value of h
STEPS       = 10000000      ##The number of steps
PRECISION   = 10**(-10)     ##the precision for Delta E / E
LAMBDA      = 0.1             ##the lambda parameter in the potential
CUTOFF      = 100           ##the cutoff for the function blowing up
POINTCUT    = 0.8           ##the amount of points to plot
B           = 0.2           ##for the linear potential

##quantum numbers
looping = False ##this option allows you to check all pairs (n,l) in the ranges below

##if no loop these are used
qn = [[12,4]] ##this is a list containing lists of [n,l] pairs to check instead
## if loop, these are used
nMin = 0
nMax = 0
lMin = 1
lMax = 2

#parameters of the potential: b, alpha_c, C, alpha h, sigma h, spin,ang.momentum
parameters = [0.782,0.857,-0.581,0.840,0.700,1,1]
#masses and energy used
constants = [0.375,0.650,1.4299235012]#
reduced_mass = 2*constants[0]*constants[1]/(constants[0] + constants[1])
#convert fm to GeV:
c_factor = 5.068

##program options:
plotting = True     ##do I want to see the plot? if false, the options below don't matter
showPlot = True     ##do I want to see the plot now?
plotSquare = False  ##do I want the square of the wavefunction?
savePlot = True     ##do I want to save the plot?

##potential options:
##1: anharmonic oscillator
##2: linear potential
##3: Harmonic oscillator
##4: interquark potential in a meson
potential = 4

##Some helpful functions

def Weinstein_Potential(x):
    if x == 0:
        return 0.
    else:
        if parameters[5] == 0:
            spin_factor = -0.75
        else:
            spin_factor = 0.25
        linear_confining_potential = parameters[0]*x  + parameters[2] 
        colour_coulomb =  -(4*parameters[1]/(x*c_factor*3.))
        colour_hyperfine = (32.*math.pi*parameters[3]*pow(parameters[4],3)*math.exp(-pow(parameters[4]*x*c_factor,2))*spin_factor)/(math.sqrt(math.pi)*9.*constants[0]*constants[1])
        return linear_confining_potential + colour_coulomb + colour_hyperfine + pow(2*reduced_mass,-1)*parameters[-1]*(parameters[-1] + 1)*pow(x*c_factor,-2.)

##The potential
def V(r,l): ##the potentials, fun to play around with
    ## uncomment the potential you want
    if potential == 1:
        return l*(l+1)/(2*r**2)+0.5*r**2+LAMBDA*r**4 ## the angular momentum barrier, the harmonic oscillator part and the anharmonic oscillator part
    elif potential == 2:
        return l*(l+1)/(2*r**2)+B*r
    elif potential == 3:
        return l*(l+1)/(2*r**2)+0.5*r**2
    elif potential == 4:
        return Weinstein_Potential(r)
    print("invalid potential number imput")
    return 0

##The f(r) function
def f(r,l,E):
    return 2*reduced_mass*(V(r,l)-E)

##the T(x) function
def T(r,l,E):
    return STEPSIZE**2*f(r,l,E)/12.0

##Finding y_(n+1) from y_n and y_(n-1):
def findNextY(y1,y2,T1,T2,T3):
    ## here are what the variables are:
    ## y1 = u(r)
    ## y2 = u(r-h)
    ## T1 = T(r-h)
    ## T2 = T(r)
    ## T3 = T(r+h)
    ## The reason for this is so that I can calculate the values of
    ## y and T once since for the next iteration, y1 becomes y2
    ## and T3 becomes T2 and T2 becomes T1
    ## The parameters are calculated by the main program
    return ((2+10*T2)*y1-(1-T1)*y2)/(1-T3)

##Integrate the set of points using trapezoid method for a set of points
def integrate(points):
    val = points[0]**2 + points[int(len(points)*POINTCUT)-1]**2 ##adds first and last point
    for i in range(1, int(len(points)*POINTCUT)-1): ##these limits for the loop are so that the ploted section is normalized to 1
        val += 2*points[i]**2
    val *= (STEPSIZE)*0.5
    return val
        
##Find the set of points for a given value of E and L
def findPoints(E,l,INTEGRATE=False):
    y1 = STEPSIZE**(l+1) ## u(h) = h^(l+1)
    y2 = 0 ## u(0) = 9
    points = [y2,y1] ##list to hold the points
    nextY = -1 ## placeholder for u(r+h)
    currentR = STEPSIZE ##the starting value of r
    T1 = 0   ## T(r-h), I use zero since it is multiplied by 0 initially and this prevents
             ## a division by zero error in V(r)
    T2 = T(STEPSIZE,l,E)   ## T(r)
    T3 = T(2*STEPSIZE,l,E) ## T(r+h)
    counter = 0 ## keeps track of the number of nodes
    current = int(points[1]/abs(points[1])) ## to keep track of the number of nodes
    for i in range(2,STEPS): ##for loop
        nextY = findNextY(y1,y2,T1,T2,T3) ##find the next y
        if abs(nextY) > CUTOFF:
            break
        if nextY == 0: ## the unlikely case where its exactly 0
            counter +=1
            current *= -1
        elif int(nextY/abs(nextY)) == -current: ## if the sign changes, it has to pass by 0
            counter +=1
            current *= -1
        points.append(nextY) ## add the next value of y to the points
        currentR += STEPSIZE #increment the r
        ##change the values
        y2 = y1
        y1 = nextY
        T1 = T2
        T2 = T3
        T3 = T(currentR+STEPSIZE,l,E) ## need to fine new T(r+h)
    if INTEGRATE:
        integral = integrate(points)**0.5
        return [[y/integral for y in points],counter]  ## divides all the points by the integral to normalize
    return [points,counter] ## returns the points and the number of nodes

##plotting points and saving the plot
def plot(points,n,l):
    power =1
    if plotSquare:
        power = 2
    R = [points[i]**power for i in range(1,int(len(points)*POINTCUT))] ## the 4/5 seems to work nicely to cut off the part that blows up
    r = [i*STEPSIZE for i in range(1,len(R)+1)] 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("r (natural units)",fontsize=15)
    ylabel = "wfn"
    if plotSquare:
        ylabel += " squared"
    ax.set_ylabel(ylabel,fontsize=15)
    title = "wfn"
    if plotSquare:
        title += " squared"
    title += " vs r, n = " + str(n) + ", l = " + str(l)+", potential = " + str(potential)
    plt.plot(r,R,axes= ax) ##the points
    ax.set_title(title,fontsize=20)
    plt.grid(True)
    if savePlot:
        filename = "wfn_"
        if plotSquare:
            filename += "sq_"
        filename += str(n) + "_" + str(l) + "_p" + str(potential)+ ".png"
        plt.gcf().savefig(filename)
    if showPlot:
        plt.show()

## this is the primary functions which brings everything together
## pass in the quantum number n and l
def findEigenvalue(n,l):
    nodesWanted = n
    minE = 1 ## after plotting a couple tests, this seems to be below the ground state
    currentMinN = findPoints(minE,l)[1]
    maxE = (n+2)**2 ## seems to work for n at least up to 100 (just in case)
    currentMaxN = findPoints(maxE,l)[1]
    answer = 0 ## placeholder for the energy eigenvalue
    answerErr = 0 ## placeholder for the error of the energy eigenvalue
    ## quick check
    while currentMinN > nodesWanted:
        minE *= 0.1
        currentMinN = findPoints(minE,l)[1]
    while currentMaxN < (nodesWanted+1):
        maxE *= 10
        currentMaxN = findPoints(maxE,l)[1]

    difference = maxE-minE ## this will be the difference between the next energy to check

    ## The strategy is to do a binary search with the condition that the last node be in the range
    ## of the two bounds.

    while True:
        tempE = 0.5*(minE+maxE)
        difference *= 0.5
        tempNodes = findPoints(tempE,l)[1]
        if tempNodes >= (nodesWanted+1):
            maxE = tempE
            currentMaxN = tempNodes
        elif tempNodes <= (nodesWanted):
            minE = tempE
            currentMinN = tempNodes
        if difference/tempE < PRECISION:
            answer = tempE
            answerErr = difference
            break
    return [answer,answerErr] ##return a combination of the error and the error as a list

######################################################################

##this is where the magic happens (
#this calls all the other functions with the proper quantum numbers and take care of printing and plotting
def makeTheMagicHappen(n,l):
    print("Calling the functions; n =",n,", l =", l)
    eigen = findEigenvalue(n,l)
    temp = "E_" + str(n)+str(l) ##just for printing
    print(temp,"=",eigen[0], "+-", eigen[1])
    points = findPoints(eigen[0],l,True)[0]
    if plotting:
        plot(points,n,l)

##This is where the magic gets called in the script
if looping:
    for n in range(nMin,nMax+1):
        for l in range(lMin,lMax+1):
            makeTheMagicHappen(n,l) ##making it happen
else:
    for p in qn:
        n = p[0]
        l = p[1]
        makeTheMagicHappen(n,l) ##making it happen

input("Done")

