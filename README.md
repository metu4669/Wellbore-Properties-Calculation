# Wellbore-Properties-Calculation

This report is about writing a MATLAB code for calculation of some properties of wellbore and these calculations will be done along the wellbore and we can make some assumption like temperature and pressure will increase linearly from wellhead to bottom hole. At the beginning of this homework, we were given with some initial values and these values can be determined based on our student id number. In my group, we used one of our student id number and the id number is xxxx.
Based on the student id number, we calculated our input data and it is given below;
 
After we found our input data, we will use them to find critical temperature and pressure value and to do this will use some correlations.
Firstly, we will use Sutton Correlations to find critical properties. I created a code segment which can calculate critical properties for Sutton Correlations.
 
‘SPpc’ and ‘STpc’ are the critical pressure and temperature respectively. We will use these values in Wichert and Aziz correlations to correct them.
 
In this homework, we were asked to draw;
•	Density vs Depth (Along wellbore)
•	Viscosity vs Depth (Along wellbore)
•	Pseudo Pressure vs Depth (Along wellbore)
To draw them, we need to find Z Factor but when we start to move from wellhead to bottom hole, temperature and pressure will increase and eventually we will have different Z Factor along the wellbore. Therefore, density, viscosity and pseudo pressure will change along the wellbore.
We divided the wellbore with a certain number and we identified the division number with ‘k’ in the code.
 
In this part, we created division and after that we wrote a code segment to calculate Apparent Molecular Weight and we created R number based on Rankine and psi value.
In this segment also, we created some matrix to store our wellbore temperature and pressure values.
First part contains, temperature and pressure matrix will hold the temperature and pressure values of each segment of wellbore. After that, by dividing them with pseudo critical and temperature values we get reduced temperature and pressure values and we showed them with ‘TemperatureR’ and ‘PressureR’ respectively.
To calculate the Z Factor, we used Redlich Kwong equation and we run this equation for each segment.
 
This code segment will calculate Z Factor along the wellbore. After that, we stored these values in Z Factor matrix and this matrix has the same amount of row as the segment number.





At this step, we will calculate viscosity by using Lee-Gonzales and Eakin Correlation. At this calculation step, we created 4 matrix and each matrix will store the related values used in Correlation and the calculation method be run for each segment. Therefore, we created loop for this calculation purpose.
 
Also, in the calculation of viscosity we will use density. Density is a function of temperature and pressure. Therefore, we will have different density value for each segment along the wellbore. Density can be calculated by
Density = Pressure * Apparent Molecular Weight / (Z Factor* R Value * Temperature)
We need to apply this equation for each segment. Therefore, before calculation of viscosity, we put density calculation code. This code will be run for each segment and result of this density will be used in viscosity calculation.
After we calculated the viscosity and density, we need to calculate pseudo pressure.
 
Before we started to calculate the pseudo pressure, we created some matrix which can hold temporary values of pseudo pressure variables during calculation.
•	l = segment of pseudo pressure calculation
•	P = temporary pressure
•	U = temporary viscosity
•	Z = temporary Z Factor
•	Den = temporary density
•	Del = Delta Pseudo Pressure
 
In this part, we created a loop which will calculate the pseudo pressure along the wellbore and write it in the ‘Pseudo’ matrix which was created to hold pseudo pressure along the wellbore. In this loop, Z Factor, density and viscosity were recalculated for pseudo pressure calculation. The reason is when we calculate the pseudo pressure we kept constant the temperature and make calculation by changing pressure from zero value to pressure of each segment. Therefore, each segment calculation were made by creating loop. At the final step of pseudo pressure calculation, we used trapezoidal rule in our calculation. 
 
Drawing Part
Drawing part includes;
•	Z Factor vs Segment Number
 
•	Density vs Segment Number
 
 
•	Pseudo Pressure vs Segment Number
 
•	Viscosity vs Segment Number
 

 
Comparing Part
 
In this code segment, we created a loop which will compare the Z Factor value of our and Z Factor value of PVTProps. In both Z Factor calculation, we used same segment pressure and temperature value and compared them. 
 
Green line is our Z Factor values and blue line is PVT prop values. We can see easily say that when pressure increase from wellhead to bottom-hole Z Factor values of ours are deviating from the PVT prop values. Therefore, Redlich Kwon Z Factor values are sensitive for lower pressure and temperature values. 
