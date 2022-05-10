DOWNLOADING AND COMPILING KNOSOS

Download the up-to-date version of KNOSOS from

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT

Specifically, the sources of KNOSOS can be found at

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Sources

This folder has all you need for compiling KNOSOS stand-alone (i.e., out of STELLOPT). There is a file with the "main" routine, knosos.f90, and a Makefile (both are ignored by STELLOPT)

You do not need to compile STELLOPT to run KNOSOS (although you can because, compiling STELLOPT will generate an executable xknosos).

The makefile has options for different computer clusters. If you need to edit it, you may find a few tips in the user manual

https://github.com/joseluisvelasco/KNOSOS/blob/master/MANUAL/KNOSOSManual.pdf


COMPUTING GAMMA_ALPHA WITH KNOSOS

Examples can be found at

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples

Specifically, the data at

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples/Gamma_alpha

allows to compute a profile of Gamma_alpha for the KJM configuration of W7-X with <beta=4%> and parabolic beta-profile taken from the VMEC RESTAPI

- Copy the content of

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples/Gamma_alpha/Input

in your working directory.

- Execute KNOSOS. This can e.g. be done by typing

../../Sources/knosos.x

if your are running from KNOSOS/Examples/Gamma_alpha/

- The calculation takes ~5-10 minutes (if compiled with MPI, the 10 flux-surfaces will be computed in parallel and it should take ~ 0.5-1 minute). This could probably be accelerated, if needed.

- The output can be compared with that of (large fort.* files have not been uploaded to GH)

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples/Gamma_alpha/Ouput

and can be processed with simple gnuplot scripts at

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples/Gamma_alpha/Scripts

The resulting pdf files should be compared to those of

https://github.com/PrincetonUniversity/STELLOPT/tree/CIEMAT/KNOSOS/Examples/Gamma_alpha/Figures

There is a one-to-one equivalence between *.plt files at Gamma_alpha/Scripts and *.pdf files at Gamma_alpha/Figures

    - "gammacsa.plt" reads file "gammacs.map" and produces a figure comparable to Figure 5 right of http://fusionsites.ciemat.es/jlvelasco/files/papers/velasco2021prompt.pdf (it uses data from a single flux-surface, which can be changed by editing the script).

    - "fraction_lambda.plt" reads file "prompt.lambda" and produces a figure similar to Figure 8 right of the same paper (the flux-surfaces read can be changed; selecting data with $6==0 gives you no smoothing, S6==1 , $6==2, ... $6=8 gives the same data with different leves of smoothing by means of moving average)

   - "gamma_alpha.plt" reads file "stellopt.knosos" and produces a radial profile of Gamma_alpha. "stellopt.knosos" contains other quantities, such as the fraction of trapped ions (last column)

- For your calculations, replace the file "boozer.txt" with your own equilibrium in boozer coordinates. At the moment, KNOSOS recognizes the format of the file from its name, so you will need to rename your file or create a link:
- "boozer.txt" (txt file employed basically at IPP Greifswald).
- "boozmn.nc" (netcdf output of BOOZER_XFORM).
- "boozmn.data" (binary output of old versions of BOOZER_XFORM).
Bear in mind that KNOSOS looks for a magnetic equilibria in the working folder and one folder above. On file "STDOUT" you can make sure that KNOSOS has read the equilibrium that you wanted it to and not another file that you left in the folder.

- The list of flux-surfaces to be computed can be modified at input.surfaces. However, I do not advice to go for a very high radial resolution, as the profile may get noisy.





