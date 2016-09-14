Assignment #3: Ray tracing

FULL NAME: Elizabeth Chu


MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  yes

2) Ray tracing sphere                     yes

3) Triangle Phong Shading                 yes

4) Sphere Phong Shading                   yes

5) Shadows rays                           yes

6) Still images                           yes
   
7) Extra Credit (up to 20 points)

10 pts — soft shadow
* DISABLE OR ENABLE SOFT SHADOW PROPERTY by changing line 56,
  which defines a SHADOW_OFFSET. Comment out line 56 if you don't
  want soft shadows. The SHADOW_OFFSET is a float and defines 
  the distance between each duplicated light source. Recommended
  value 0.02 ~ 0.2
* Implemented soft shadow by simulating area light, changing each 
  light source from:
                                   +     +     +
                                     +   +   +  
                                       + + +    
          +             to         + + + + + + +
                                       + + +    
                                     +   +   +  
                                   +     +     +
* Each light's color intensity is attenuated by a factor of 25. 
* Still images are provided for both the hard shadow and soft 
  shadow renderings. 


