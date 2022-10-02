# Design Goals

## For users

- Easy to switch to a new EOS model
- Easy to switch to another **Flash Calculation** algorithm
- Easy to add a new IMPES/IMPEC method

## For developers
- Easy to change to a different grid structure
- Easy to add a new spatial discretization
- Easy to add a new temporal discretization
  
## For performance
- Easy to add a new linear solution method
- Easy to parallelize linear solvers
- Easy to optimize code at certain steps

## Overall design

<img src="./OCPStructure.png" alt="Structure" width="500"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## Top Structure

<img src="./TopStructure.png" alt="Top Structure" width="600"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## Solve Structure

<img src="./SolveStructure.png" alt="Solve Structure" width="600"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## WorkFlow

<img src="./WorkFlow.png" alt="WorkFlow" width="600"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## Linear solution methods

<img src="./OCPLinearSolver.png" alt="Linear solver" width="500"
     style="display: block; margin: 0 auto"/>