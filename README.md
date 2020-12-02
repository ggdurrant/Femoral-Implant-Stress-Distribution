# Femoral Implant Stress Distribution
This project analyzes the difference in stress distributions on bone caused by implants by applying asymmetric composite beam theory. Implants like hip replacements can greatly affect stresses in the surrounding bone, so location of the implant and bone density must be accounted for to ensure a solid fit and reduce complications. Theoretical moments are applied to bone, with the moments and type of implant held constant so just the differences in bone and implant location can be examined.  

## Cases
For each of the types of bone:
 - Homogenous young bone, elastic modulus assumed to be 17 GPa
 - Homogenous old bone, elastic modulus assumed to be 13 GPa
 - Heterogenous bone, elastic modulus assumed to be a function of bone mineral density (BMD), calculated from bone scan
We then examine a constant steel implant, elastic modulus of 200 GPa at 4 locations:
 - a
 - b
 - c
 - d
 
![Locations](https://github.com/ggdurrant/Femoral-Implant-Stress-Distribution/blob/main/images/implantlocs.PNG)

The implant and BMD weighted centroids are calculated in each case, such as:
![Centroid](https://github.com/ggdurrant/Femoral-Implant-Stress-Distribution/blob/main/images/centroid.png)

## Results
The max stresses and corresponding locations for each case can be calculated, as well as looking at the distribution of compressive and tensile stresses, which are plotted with the neutral axis. 

### The stress distribution of the variable BMD heterogenous bone 
### without an implant:
![NA](https://github.com/ggdurrant/Femoral-Implant-Stress-Distribution/blob/main/images/noImplant.png)

### with implant at location d:
![NAImplant](https://github.com/ggdurrant/Femoral-Implant-Stress-Distribution/blob/main/images/implantD.png)
