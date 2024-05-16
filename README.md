# Stress Analysis

## Description

This repository contains two files, `circapp` and `ellipseappl`, which have been created for the purpose of implementing pure tension and pure bending applications.

Initially, the active degrees of freedom are retained. Each node has two degrees of freedom, and at the left end, all nodes except one have only one active degree of freedom, while the other is inactive. In the exceptional node, both degrees of freedom are inactive, meaning that neither of the two accepts motion. For the application of force and stress, a different method will be used.

For the case of pure bending, for each element of the surface of the right end, I distribute the forces to produce the desired moment. Essentially, I use the relationship, moment / distance. Then, I distribute the above force equally into 2 parts and assign them to the 2 degrees of freedom (i.e., to the degrees of freedom of the nodes) of the surface of that specific element.

For the case of pure tension, we follow the same logic, as we want to distribute the force equally to all elements, and then for each surface of each element corresponding to a part of the right end of the specimen, we distribute this force equally to the middle of the 2 nodes that are part of this surface. As a result, the two ends of the specimen will have nodes with half the force applied compared to the force applied to the intermediate nodes.




