all: input
	fluidity  -v5 balzano1_2plus1.flml &> output.txt

input: clean
	cd mesh; ./makemesh.sh

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.png
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* 
config:
	diamond -s ../../tools/fluidity_options.rng  balzano1_2plus1.flml 

