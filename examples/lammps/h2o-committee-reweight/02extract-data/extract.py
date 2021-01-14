# collect of misc function
import argparse
import sys
import json
import numpy as np
import matplotlib.pyplot as plt


def ipi_extra(arg):
    """
    function to extract the scaled potential from ipi extra file
    """
    extra_file = arg.ipi_extra_file
    print("find file {0}".format(extra_file))
    pot_num = arg.pot_num
    print("we have {0} potentials".format(pot_num))
    if arg.factor :
        unit_conversion = arg.factor
    else:
        unit_conversion = 1.0

    print("extracting the potential")
    fr = open(extra_file, "r")
    frame_pot_list = []
    for idx, line in enumerate(fr.readlines()):
        if (idx%2 == 0):
            continue
        else:
            info = json.loads(line)
            committee = info["committee"]
            ml_pot = []
            for i in range(pot_num):
                ml_pot.append(committee[i]["v"])
            frame_pot_list.append(ml_pot)
    fr.close()
    #save file
    frame_pot_list = np.array(frame_pot_list)*unit_conversion
    np.savetxt("ml_pot.dat", frame_pot_list)
    print("save file to ml_pot.dat")
    print("FINISHED")

def main():
    print("Description\n------------")
    parser = argparse.ArgumentParser(description="""
    tools to reformat the i-pi .extra data""")

    #subparsers = parser.add_subparsers()

    # tools
    parser.add_argument("ipi_extra_file", type=str,
            help="ipi raw extra file")
    parser.add_argument("pot_num", type=int,
            help="the number of potential in committee")
    parser.add_argument("--factor", nargs="?", const=1.0,
            type=float, help="conversion factor used to change the unit of energy")
    parser.set_defaults(func=ipi_extra)



    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)

    args.func(args)

if __name__ == "__main__":
    main()
