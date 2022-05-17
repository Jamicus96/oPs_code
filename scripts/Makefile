RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: Hist_Time_Profile.exe

Hist_Time_Profile.exe: Hist_Time_Profile.cpp
	g++ -g -std=c++1y Hist_Time_Profile.cpp -o Hist_Time_Profile.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
