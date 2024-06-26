# Preliminary-OD

MAE 180A Orbit Determination Code Instructions

Created by: Miles Puchner

3/17/23



CONTENTS:
--> "READ ME" .txt file containg code description and instructions

--> "MAE180MasterLiveScript" .mlx file to quickly run code using the given data

--> "OrbitComp" .m file contianing the main orbital determination code

--> "figures" .m file containing the orbital plotting code

--> "MAE 180A Final Project Writeup" .pdf file containing the report



INSTRUCTIONS:
--> Open the "MAE180MasterLiveScript" file to run the orbital determination method

--> Within this file all input variables are already intialized according to the input data given to our group

--> The 'model' varaible may be tweaked according to the instructions on the "MAE180MasterLiveScript" file in 
    order to run the main "OrbitComp" function using either an Oblate-Spheroid Earth model or a Spherical Earth model
   --> The 'model' variable is set to 'oblate' by default

--> Run the "MAE180MasterLiveScript" file to determine the initial and propagated states of the given satellite data
    and to plot these outputs

--> If MATLAB produces an error, ensure the specified pathway within MATLAB leads to the .zip file containing this 
    "READ ME" .txt file

--> The "OrbitComp" and "figures" functions will automatically run once the "MAE180MasterLiveScript" is ran

--> For specifics on the orbital determination and plotting methods used, open the respective
    .m function files contained with the .zip file
   --> At the beginning of each function ("OrbitComp" or "figures") is an in-depth description
   --> The functions are also well commented for ease of understanding

--> If both the spherical and oblate-spheroid Earth model data is desired, the "MAE180MasterLiveScript" file should
    be ran twice with the 'model' variable changed between 'oblate' and 'sphere'.



FUNCTION DESCRIPTIONS:
--> The "OrbitComp" function utilizes Gauss' Angles Only method of preliminary orbital determination to define the 
    Keplerian orbital elements and state vectors at an intial and propagated state. It can handle both a spherical Earth
    model as well as an oblate-spheroid Earth model which is specified by the input variable 'model'.
   --> NOTE: The "OrbitComp" function employs an iterative refinement method when calculating slant values. By default, 
       the loop will terminate when one of the following two conditions is met: the newly calculated slant values rounded to
       ten decimal places equals the previous iteration's slant values rounded to ten decimal places, OR, one hundred
       iterations have been achieved.

--> The "figures" function utilizes a standard ECI to orbital elements algorithm to create the following two plots. 
   --> Given the initial and propagated state vector of the main satellite and the type of Earth model used in their 
       calculations, the function will output a plot of the initial and propagated orbits and positions at the respective 
       epochs defined by each state vector. If the model specified is the Spherical Earth model, this plot will only contain
       one orbit as the orbit does not change in 3 dimensional space for this model.
   --> Given the propagated state vector of the main satellite and the state vector of another satellite corresponding to the
       same epoch, the function will output a plot of the propagated main satellite's orbit and position as well as the 
       orbit and position of the compare satellite.
