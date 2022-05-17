#######################################################################################
# Functions demonstrating the use of the TrackNav, TrackCursor and TrackNode classes
# All require a ROOT file generated with tracking information ON.
#
# see https:##github.com#snoplus#rat#blob#master#example#macros#tracking#tracking.mac
#     http:##snopl.us#docs#rat#user_manual#html#node214.html
#
# Particle tracks are divided into a series of backward looking 'Nodes' containing
# the step information. The TrackCursor class controls the navigation through these
# nodes and between tracks.
#
# Secondary particles appear new 'child' tracks at the node that created them.
#
# Primary particles are the children of an inital track with a single node
# representing the event
#
#
# J dunger - 2014-08-20 <jack.dunger@ox.ac.uk> : New file
#######################################################################################

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import ROOT
import rat
import argparse

def set_parser():
    # Define arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='Input args')

    parser.add_argument('--sim_file', '-f', type=str, default=None,
                        dest='sim_file', required=True,
                        help='Root file from simulation')
    parser.add_argument('--which_func', '-w', type=int, default=None,
                        choices=[0, 1, 2],
                        dest='which_func', required=True,
                        help='Which function to run')

    return parser.parse_args()

def main():
    args = set_parser()

    functions = [plot_primary_mom, plot_primary1_disp, find_secondaries]

    # Run desired function
    functions[args.which_func](args.sim_file)


def plot_primary_mom(filename):
    """ Plot the primary particle momenta for a set of events.

    :param file_name: Path to the RAT DS file
    """
    hist = ROOT.TH1D("", "", 10, 0, 5)

    # loop over events
    for ds, _ in rat.dsreader(filename):
        nav = ROOT.RAT.TrackNav(ds)
        cursor = nav.Cursor() # cursor at head node of track
        # head node is a track of its own
        # real tracks for generated particles are its children

        for child in range(0, cursor.ChildCount()):
            node = cursor.Child(child) # grab node at start of child track
            hist.Fill(node.GetMomentum().Mag()) # nodes contain particle information
            print(node.GetParticleName())

    return hist

def plot_primary1_disp(filename):
    """ Plot the 1st primary particle's displacement.

    :param file_name: Path to the RAT DS file
    """
    hist = ROOT.TH1D("", "", 10, 0, 1000)

    # event loop
    for ds, _ in rat.dsreader(filename):
        nav = ROOT.RAT.TrackNav(ds)
        cursor = nav.Cursor()

        if cursor.ChildCount():
            cursor.GoChild(2) # move cursor to start of child track
            node = cursor.Here() # grab the node here
            init_pos = node.GetPosition()

            node = cursor.GoTrackEnd() # returns correct node AND moves cursor
            fin_pos = node.GetPosition()

            disp = (fin_pos - init_pos).Mag()
            hist.Fill(disp)

    return hist

def find_secondaries(filename):
    """ Find the processes with secondaries.

    :param file_name: Path to the RAT DS file
    """
    # event loop
    for ds, _ in rat.dsreader(filename):
        nav = ROOT.RAT.TrackNav(ds)
        cursor = nav.Cursor()

        # primary particle loop
        for child in range(0, cursor.ChildCount()):
            cursor.GoChild(child)

            # # Step along particle track
            # while not cursor.IsTrackEnd():
            #     node = cursor.Here()
            #     if cursor.ChildCount():
            #         print(node.GetProcess())
            #     cursor.GoNext()

            cursor.GoTrackEnd()
            node = cursor.Here()
            print(node.GetProcess())

            cursor.GoTrackStart()
            cursor.GoParent()


if __name__ == '__main__':
    main()