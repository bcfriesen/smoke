CXX=CC -O3 -ipo -g -qopt-report=5 -qopt-report-phase=vec,cg -std=c++11

TARGET = smoke

OBJECTS = main.o grid.o transport.o spectrum.o radioactive.o rng.o

$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS)

.cc.o: 
	$(CXX) -c -o $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS) *.oo

realclean: clean
	rm -f $(TARGET)
