# mcEvol

The installation should be straighforward, the only external dependence is on ROOT which is used for filling histograms with pdfs or kernel.
1) Clone the repository
`git clone git@github.com:zleba/mcEvol.git`
2) Please check that the paths to the ROOT are set correctly in the Makefile
3) Please check the compiler setting (variable CC) in the Makefile. At least g++ 4.7 is required.

For running simply call `./mcEvol` which will fill the histograms (in root file) with z and PT spepctrum at various scales.
The "event loop" is in the main fuction in src/main.cpp, especially one can set here the number of "events"

During short init period the splitting functions are transormed into splines, therefore the running at LO, NLO or NNLO takes roughly the same time.

Generation of 1 million "events" takes about 20s, including the initialisation.
