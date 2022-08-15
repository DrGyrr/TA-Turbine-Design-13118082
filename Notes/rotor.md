# Rotor Module Note

## Rotor Computing Funtions

### ComputeR1: (deprecated)
This function is deprecated, thus not able to be used anymore. But this section will serve as documentation nevertheless

=>INPUT: 10x flow coeff, 10x work coeff, cycle number  
> \* flow coeff = Cm4/U4

=>OUTPUT: Effts


### ComputeR2:
This function takes definition of flow coeff as Cm5/U4  
> \* flow coeff = Cm5/U4

INPUT: 
- 10x flow coeff, 
- 10x work coeff, 
- cycle number,
- gparamset number
- rpm number
- - [ ] to do
 

>\* loop => Cm4.<br>
This function use loop to compute Cm4

OUTPUT: Dictionaries of categorized computed variables

### ComputeR3 
This function takes definition of flow coeff as Cm4/U4  
> \* flow coeff = Cm4/U4  

INPUT: 
- 10x flow coeff, 
- 10x work coeff, 
- cycle number,
- gparamset number
- rpm number

OUTPUT: Dictionaries of categorized computed variables
>\* loop => Cm4.<br>
This function use loop to compute Cm4


### ComputeR4