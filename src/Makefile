
#CXX=CC -O3 -fastsse -static

TARGET = smoke

OBJECTS = main.o grid.o transport.o spectrum.o radioactive.o

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
