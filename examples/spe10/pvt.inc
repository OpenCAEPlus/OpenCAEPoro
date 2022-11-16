PVDO      -- pvt for live oil
-- Pressure     Bo      viscosity 
      300      1.05       2.85      
      800      1.02       2.99      
     8000      1.01       3.00          
/                            

-- Rs remains constant after Pbub = 3824.32 psi at the PVTO table 



PVTW
    6000    1.01    3.00E-06  0.3	0 
/

     
ROCK       -- reference Pressure and rock compressibility
LINEAR01          6000 1E-6  
/

DENSITY   -- oil wat gas @surface(lbm/scf)
        53.0  64.00   0.01078
/
