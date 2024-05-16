# Stress Analysis

## Description

This repository contains files `circapp` and `ellipseappl` designed for implementing pure tension and pure bending applications, respectively.

### Tension Implementation:

In the case of tension, forces are evenly distributed across the elements. Each surface of an element at the right end of the specimen receives an equal share of the force, which is then equally distributed to the two nodes forming that surface. As a result, nodes at the ends of the specimen experience half the force compared to intermediate nodes.

![Tension Implementation](https://github.com/Vasilisdi/stress-analysis/assets/24864439/78fa486f-0074-4224-82d2-28f171b03579)

### Bending Implementation:

For pure bending, forces are distributed to generate the desired moment on each element's surface at the right end of the specimen. The force is equally divided into two parts, assigned to the degrees of freedom of the nodes forming that surface.

![Bending Implementation](https://github.com/Vasilisdi/stress-analysis/assets/24864439/be1e1195-187c-467e-8be2-dc4c422c7ecb)

### Matrices and Vectors:

- Single Element Surface Matrix:
  ![Surface Matrix](https://github.com/Vasilisdi/stress-analysis/assets/24864439/539c355c-1bb0-4b16-b688-a302a6ab8fb8)

- Strain Matrix:
  ![Strain Matrix](https://github.com/Vasilisdi/stress-analysis/assets/24864439/56c96ba8-d18c-48ff-a0d7-580748b1dd23)

- Stress Vector:
  ![Stress Vector](https://github.com/Vasilisdi/stress-analysis/assets/24864439/6b314317-fb9c-4ff3-bf62-aaf1469d0a82)

### Stress Results:

- Tension Stress Results:
  ![Tension Stress](https://github.com/Vasilisdi/stress-analysis/assets/24864439/30add36f-6de3-465b-8683-eb2e7caa844e)

- Bending Stress Results:
  ![Bending Stress](https://github.com/Vasilisdi/stress-analysis/assets/24864439/30add36f-6de3-465b-8683-eb2e7caa844e)


