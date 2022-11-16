RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: Hist_Time_Profile.exe getPDF.exe ClassifierResults.exe

Hist_Time_Profile.exe: Hist_Time_Profile.cpp
	g++ -g -std=c++1y Hist_Time_Profile.cpp -o Hist_Time_Profile.exe ${RAT_CONFIG}

getPDF.exe: getPDF.cpp
	g++ -g -std=c++1y getPDF.cpp -o getPDF.exe ${RAT_CONFIG}

ClassifierResults.exe: ClassifierResults.cpp
	g++ -g -std=c++1y ClassifierResults.cpp -o ClassifierResults.exe ${RAT_CONFIG}

get_evt_info.exe: get_evt_info.cpp
	g++ -g -std=c++1y get_evt_info.cpp -o get_evt_info.exe ${RAT_CONFIG}
	
clean :
	rm *.exe
