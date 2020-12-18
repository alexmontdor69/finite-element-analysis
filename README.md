## FEA Solver -

Finite Element Analysis

4 basic types of Element

- Truss Element
- beam
- place 3 Nodes
- plane 4 nodes

## Data File

- Define the number of Nodes
  <code>NODE,4</code>

- Define the number of Elements
  <code>ELEM,3</code>

- Define the coordinates of Node : Node Id, position X, position Y
  <code>NBLOCK
  1,0,1000
  2,2000,0
  3,1000,0
  4,00,00
  ENDNODEBLOCK</code>

- Define the element
  <code>EBLOCK</code>

- ET give the type of element
  <code>ET,1</code>

- Descrition of the element : element Id, material Id, cross area, Node Id A, Node Id B
  <code>3,1,1000,1,2</code>

<code>ET,2</code>

- Descrition of the element : elementId, material, cross area, second moment of area / inertia, node, node
  <code>1,1,5000,1000000,3,2
  2,1,5000,1000000,4,3
  ENDEBLOCK</code>

- Material Description : K for Young Modulus; v for poisson coefficient
  <code>MAT
  K,207000
  EMAT</code>

- Forces description : #node, fx, fy, momentz
  <code>BForces
  3,1000,-2000,0
  ENDBForces</code>

- Displacement description : #node, u, v, rotz ; 0 = blocked, 1=available
  <code>BDis
  1,0,0
  4,0,0,0
  ENDBDis</code>

## command line

### compilation

g++ -o sofea1 sofea.cpp
or
gcc sofea.cpp

Info
gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04)
