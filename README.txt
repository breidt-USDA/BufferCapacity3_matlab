BufferCapacity3 is a Matlab program developed in Matlab version 2022a/b

For many acidic foods, including fermented and acidified vegetables, salsas, salad dressings, and others, maintaining a pH below 4.6 is a critical control to prevent botulism. The pH of acidic foods is controlled by acid content, low acid ingredients, and buffering; however, buffering of foods remains largely uncharacterized. A Matlab GUI program BufferCapacity3 was developed to automate the process of quantifying the buffers present in foods utilizing acid/base titration data. The BufferCapacity3 program may be used to aid product development and help assure pH control and safety of acidic foods.

The app may be installed in Matlab by drag-and-drop of the BufferCapacity3.mlappinstall into the Matlab Command Window.

The BufferCapacity3.mlapp file contains the Matlab AppDesigner code to build the executable using that program.  

Example input files are included with data from a Hanna Instruments titrator (Model 931, Hanna Instruments, Smithfield, RI, USA). The data consist of two titration files for an acid (HCl) and base (NaOH) titration of a cucumber juice medium:

Sample input files included:

EXAMPLE_ACID.RPT 
EXAMPLE_BASE.RPT

Sample output files included:

Ingredient.mat: the ABD class structure containing most program output
Ingredient_BCdatapoints.csv: a comma delimited spreadsheet file with the original buffer capacity curve
Ingredient_BCmodelcurve.csv: a comma delimited spreadsheet file with the BC model buffer curve
Ingredient_buffertable.csv: A comma delimited spreadsheet file with the concentraton, pK data for the buffer table
Ingredient_graph.pdf: a pdf file showing the buffer model graph
Ingredient_info.csv: a comma delimited spreadsheet with summary details of the buffer model, including the ingredient name, concentation, titration results, etc.
Ingredient_paramfile.csv: a copy of the parameter file used for generating titration data

Additional files include:
 
paramfile.csv: The parameter file used by the program for processing titration files
help.pdf: a pdf help file with a link to software information

The dependency folder includes all of the Matlab functions needed for the AppDesigner code and must be in the Matlab path when using AppDesigner to run the program. 

USDA/ARS work is not subject to copyright, and is supplied with the Creative Commons Zero (CCO) license.

Modeling references:

Breidt F, Skinner CR. 2022. Buffer models for pH and acid changes occurring in cucumber juice fermented with Lactiplantibacillus pentosus and Leuconostoc mesenteroides. J Food Prot 85(9):1273-1281. https://doi.org/10.4315/JFP-22-068.

Price RE, Longtin M, Conley Payton S, Osborne JA, Johanningsmeier SD, Bitzer D, Breidt F. 2020. Modeling buffer capacity and pH in acid and acidified foods. J Food Sci 85(4):918-925. https://doi.org/10.1111/1750-3841.15091.

Longtin M, Price RE, Mishra R, Breidt F. 2020. Modeling the buffer capacity of ingredients in salad dressing products. J Food Sci 85(4):910-917. https://doi.org/10.1111/1750-3841.15018.

Author contact: 
Fred Breidt, PhD
USDA/ARS Research Microbiologist
322 Schaub Hall, Box 7624
NC State University
Raleigh, NC 27695-7624
Tel: 919-513-0186
https://www.ars.usda.gov/southeast-area/raleigh-nc/fsmqhru/
fred.breidt@usda.gov

   

