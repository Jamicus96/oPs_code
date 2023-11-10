RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux

all: get_evt_info.exe make_hists.exe compare_hists.exe printClass_ntuple.exe

get_evt_info.exe: get_evt_info.cpp
	g++ -g -std=c++1y get_evt_info.cpp -o get_evt_info.exe ${RAT_CONFIG}

make_hists.exe: make_hists.cpp
	g++ -g -std=c++1y make_hists.cpp -o make_hists.exe ${RAT_CONFIG}

compare_hists.exe: compare_hists.cpp
	g++ -g -std=c++1y compare_hists.cpp -o compare_hists.exe ${RAT_CONFIG}

printClass_ntuple.exe: printClass_ntuple.cpp
	g++ -g -std=c++1y printClass_ntuple.cpp -o printClass_ntuple.exe ${RAT_CONFIG}
	
clean :
	rm *.exe