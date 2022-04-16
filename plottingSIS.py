
# ************** SIS Model differential equation Solver vs. 1.01 May 20, 2020
# miscellaneious constants and stuff
Title  = "SIS Model Differential Equation Solver"
Vers   = "vs. 1.01 May 20, 2020  David A. Larrabee"
EOL = "\n"
Header = " day, Susceptible, infectious,   formula,     error,    % infected "
e=2.7182818284  # Euler's number (approx)

# ***************************************************************************
#   Imports and definitions
# **************************************************************************
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def f(y, t, parm):             # returns the derrivative of f(t) (y)
    I = y[0]                   # get the current value of y
    plam,gammar, Npop = parm             # get the parameter
    derivative = [plam*I*(1-I/Npop)] -gammar*I     # calculate the derrivative
    return derivative

#****************************************************************************
#   These are parameters the User can Input
#****************************************************************************
plam = 0.2587      # growth rate at start (infections per day per infector)
gammar = 0.045     # the daily probability of recovery from the virus
Npop = 58500       # Total "population" (Hubei in thousands)
I0 = 1.619         # Initial number of infected (hubei in thousdands)
Finaltime = 90.    # The simulation goes from time = 0 to here (in days)
Increment = 1.0    # How often you want the result retured (in days)
#****************************************************************************
#     Set up for the simulation
#****************************************************************************
print(Title)
print(Vers)
print(" ")

arrayS = [0 for a in range(90)]
parm = [plam,gammar, Npop]   #ODE solver wants a list of Parameters
y0 = [I0]       #ODE also wants a list of initial values

# Make an array of the times at which we want a solution
times = np.arange(0.0, Finaltime, Increment)

#*************************************************************
#             Call the ODE solver
#*************************************************************
Solutions = odeint(f, y0, times, args=(parm,))

#*******************************************************************
# Compare the results to theory and print out the results to screen
#*******************************************************************

print(Header)
lambdastar = plam - gammar
Nstar=Npop*lambdastar/plam

for i in range(len(times)): # printout and compare with actual solution
    t = times[i]
    I = Solutions[i,0]
    S = Npop - I
    arrayS[i] = S
    val = Nstar/(1+(Nstar-I0)/I0 * e**(-1.0*lambdastar*t) ) # analytic solution
    error = (val-I)/val
    pinfect = I/Npop*100.0
    print("{:4.1f}".format(t), "{:11.1f}".format(S),"{:11.1f}".format(I),
          "{:11.1f}".format(val),"  ","{:+8.4E}".format(error),
          "  ","{:8.4f}".format(pinfect) )
    plt.plot(times, arrayS)
    plt.plot(times, Solutions)


print("All Done")
plt.title("SIS")
plt.xlabel("Giorni")
plt.ylabel("Persone")
plt.show()
