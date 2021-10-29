(TeX-add-style-hook "KNOSOSManual"
 (lambda ()
    (LaTeX-add-bibliographies
     "./bibliography")
    (TeX-add-symbols
     '("ivar" 7)
     '("ivarl" 7)
     '("vlink" 1)
     '("mylink" 1)
     '("mytarget" 1)
     '("makesub" 1)
     '("fsa" 1)
     '("tocheck" 1)
     '("todo" 1)
     '("vec" 1)
     '("matr" 1)
     "notf"
     "true"
     "false"
     "stella"
     "DKES"
     "JMP"
     "VMEC"
     "SFINCS"
     "EUTERPE"
     "FORTEC"
     "BOOZERXFORM"
     "NEO"
     "MPI"
     "FFTW"
     "PETSC"
     "bb"
     "bB"
     "bF"
     "bj"
     "bk"
     "bPi"
     "bu"
     "bE"
     "es"
     "epol"
     "etor"
     "epolt"
     "cL"
     "vpar"
     "bun"
     "bv"
     "bR"
     "dd"
     "bnabla")
    (TeX-run-style-hooks
     "xstring"
     "hhtensor"
     "tcolorbox"
     "float"
     "bm"
     "hyperref"
     "url"
     "fancyhdr"
     "natbib"
     "square"
     "color"
     "usenames"
     "verbatim"
     "longtable"
     "tabularx"
     "graphicx"
     "amsmath"
     "amssymb"
     "latex2e"
     "rep11"
     "report"
     "11pt"
     "twoside"
     "singlespace"
     "cover"
     "preface"
     "introduction"
     "installation_and_run"
     "input"
     "output"
     "diagnostics")))

