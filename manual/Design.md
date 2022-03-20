# Design Goals

## For users

- Easy to switch to a new EOS model
- Easy to switch to another **Flash Calculation** algorithm
- Easy to add a new IMPES/IMPEC method

## For developers
- Easy to change to a different grid structure
- Easy to add a new spatial discretization
- Easy to add a new temporal discretization
  
## For optimzation
- Easy to add a new linear solution method
- Easy to parallelize linear solvers
- Easy to optimize code at certain steps

## Overall design

<img src="OCPStructure.png" alt="Structure" width="500"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## Flow chart

<img src="FlowChart.png" alt="Flow chart" width="600"
     style="display: block; margin: 0 auto"/>

<div STYLE="page-break-after: always;"></div>

## Linear solution methods

<img src="OCPLinearSolver.png" alt="Linear solver" width="500"
     style="display: block; margin: 0 auto"/>