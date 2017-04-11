import sys
# Define a parameter optimization problem in relation to the 
# LAS algorithm.
from las_declaration import LAS

from opal import ModelStructure, ModelData, Model

from opal.Solvers import NOMAD

# Return the error measure.
def total_time(parameters, measures):
    with open("measures", "a") as f:
        f.write("parameters = %r, measures = %r\n" % (parameters, measures))
    # Now with numerical integration, we want to minimize total time, not the time/rel
    # return sum(map(float, measures["TIME"])) / sum(map(float, measures["RELATIONS"]))
    return sum(map(float, measures["SIEVETIME"]))

# Define parameter optimization problem.
data = ModelData(LAS)
struct = ModelStructure(objective=total_time)  # Unconstrained
model = Model(modelData=data, modelStructure=struct)

# Solve parameter optimization problem.
NOMAD.set_parameter(name='DISPLAY_STATS',
                    value='%3dBBE  %7.1eSOL  %8.3eOBJ  %5.2fTIME %5.2SIEVETIME')
# the following limits the number of evaluations to 100
NOMAD.set_parameter(name='MAX_BB_EVAL', value=100)
NOMAD.solve(blackbox=model)
