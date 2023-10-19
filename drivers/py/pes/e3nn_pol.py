"""Interface with e3nn to run machine learning electric dipole"""

import sys
import os
# import json5 as json
import json
import numpy as np
from ase.io import read 
from py.pes.dummy import Dummy_driver
# import importlib
from miscellaneous.elia.nn.functions import get_model
from miscellaneous.elia.functions import str_to_bool
from miscellaneous.elia.functions import add_default

# "args":["-m","e3nn_pol","-a","localhost","-p","0000","-u","-o",
#            "/home/stoccoel/google-personal/miscellaneous/miscellaneous/elia/nn/water/instructions.json,/home/stoccoel/google-personal/miscellaneous/miscellaneous/elia/nn/water/LiNbO3/MD/start.xyz"],

# def str_to_bool(s):
#     s = s.lower()  # Convert the string to lowercase for case-insensitive matching
#     if s in ("1", "true", "yes", "on"):
#         return True
#     elif s in ("0", "false", "no", "off"):
#         return False
#     else:
#         raise ValueError(f"Invalid boolean string: {s}")


# def get_class(module_name, class_name):
#     try:
#         # Import the module dynamically
#         module = importlib.import_module(module_name)
        
#         # Get the class from the module
#         class_obj = getattr(module, class_name)
        
#         # Create an instance of the class
#         #instance = class_obj()
        
#         return class_obj
    
#     except ImportError:
#         raise ValueError(f"Module '{module_name}' not found.")
#     except AttributeError:
#         raise ValueError(f"Class '{class_name}' not found in module '{module_name}'.")
#     except Exception as e:
#         raise ValueError(f"An error occurred: {e}")
#     # return None

# # Usage example
# module_name = "my_module"  # Replace with your module's name
# class_name = "MyClass"     # Replace with the class name you want to allocate

# instance = create_instance(module_name, class_name)
# if instance:
#     instance.some_method()  # Call a method on the allocated instance


class e3nn_pol(Dummy_driver):

    opts_default = {
        "compute-BEC" : True ,
        "print" : True,
        "delete-files" : True,
        "save" : {
            "dipole" : "dipole.txt",
            "BEC"    : "bec.txt"
        }
    }

    def __init__(self, args=None):
        self.error_msg = """The parameters of 'e3nn_pol' are not correctly formatted. \
            They should be two strings, separated by a comma."""
        super().__init__(args)

        # self.dipole_file = "dipole.txt"
        # self.bec_file    = "bec.txt"

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in librascal
        """
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        if len(arglist) >= 2:
            # self.model_file = arglist[0] # file with the torch.nn.Module
            info_file       = arglist[0] # json file to properly allocate a 'model' object
            parameters_file = arglist[1] # *.pth file with the model parameters
            try :
                opts_file       = arglist[2] # json file with some parameters for this class
            except :
                print("\tNo options file provided: using the default values")
                opts_file = None
        else:
            sys.exit(self.error_msg) # to be modified

        print("\n\tThe driver is 'e3nn_pol'")
        print("\n\tLoading model ...")
        self.model = get_model(info_file,parameters_file)

        if opts_file is not None :

            try : 
                # Open and read the JSON file
                with open(opts_file, "r") as json_file:
                    # Parse the JSON data and save it to a variable
                    self.opts = json.load(json_file)
                self.opts = add_default(self.opts,self.opts_default)
            except :
                print("\tError reading the options file '{:s}'".format(opts_file))
                print("\tThe default values will be used")
                self.opts = add_default(None,self.opts_default)
        else :
            self.opts = add_default(None,self.opts_default)
        
        print("\tComputing BECs: {:s}".format("yes" if self.opts["compute-BEC"] else "no"))
        print("\tInitialization completed\n")

        if self.opts["delete-files"]: 
            for file in [self.opts["save"]["dipole"],self.opts["save"]["BEC"]]:
                if file is not None and os.path.exists(file):
                    print("\tRemoving file '{:s}'".format(file))
                    os.remove(file)
        
            print("\n")

        self.count = 0 

        pass
 
    def __call__(self, cell, pos):
        """Get energies, forces, stresses and extra quantities"""

        self.count += 1

        if self.opts["print"] : 
            print(" @calling 'e3nn_pol' for the {:d}th time".format(self.count))

        # Check that if 'cell' has some np.inf values (isolated system)
        # the all the other elements are zero
        has_inf_values = np.any(np.isinf(cell))
        if has_inf_values:
            # print("The array contains inf values.")

            non_inf_mask = np.logical_not(np.isinf(cell))
    
            # Check if all non-inf values are zero
            all_non_inf_are_zero = np.all(cell[non_inf_mask] == 0)
            
            if not all_non_inf_are_zero:
                raise ValueError("Error with 'cell': the the are both inf and non-zero values.\nIs this the cell of an isolated system?")
    
        # For isolated systems the diagonal elements of the cell are np.inf
        # we need to replace them to be sure that 'vir' will be zero
        cell[np.isinf(cell)] = 0.0

        # Get vanishing pot, forces and vir
        pot, force, vir, extras = super().__call__(cell,pos)
        extras = {}

        if self.opts["compute-BEC"] : 
            dipole, bec, X = self.model.get_value_and_jac(cell=cell,pos=pos)
            extras["BEC"] = bec.tolist()
            print("sum rule:",bec.sum(dim=0))
        else :
            dipole, X = self.model.get(cell=cell,pos=pos)

        extras["dipole"] = dipole.tolist()

        # # get 'dipole' and 'BEC' tensors
        # #with torch.no_grad():
        # dipole,X = self.model.get(cell=cell,pos=pos,what="dipole",detach=not self.opts["compute-BEC"])
        # dipole = dipole[0] # remove the 'batch_size' axis
        # extras["dipole"] = dipole.tolist()

        # # file = self.opts["save"]["dipole"]
        # # if file is not None and str_to_bool(file):
        # #     with open(file, 'a') as f:
        # #         line = "{:d} : {:s}\n".format(self.count,dipole.detach().numpy().tostring())
        # #         f.write(line)
        # #         #np.savetxt(f, dipole.detach().numpy(),delimiter=" ")        

        # if self.opts["compute-BEC"] :

        #     N = len(X.pos.flatten())
        #     bec = np.full((N,3),np.nan)

        #     dipole[0].backward(retain_graph=True)#.flatten() # -> row 1 of BEC.txt
        #     bec[:,0] = X.pos.grad.flatten().detach().detach().numpy()
        #     X.pos.grad.data.zero_()

        #     dipole[1].backward(retain_graph=True)#.flatten() # -> row 2 of BEC.txt
        #     bec[:,1] = X.pos.grad.flatten().detach().detach().numpy()
        #     X.pos.grad.data.zero_()

        #     dipole[2].backward(retain_graph=True)#.flatten() # -> row 3 of BEC.txt
        #     bec[:,2] = X.pos.grad.flatten().detach().detach().numpy()
        #     X.pos.grad.data.zero_()

        #     # Axis of bec :
        #     #   1st: atoms index (0,1,2...)
        #     #   2nd: atom coordinate (x,y,z)
        #     #   3rd: dipole direction (x,y,z)
        #     # bec = bec.T.reshape((-1,3,3)) 
        #     bec = bec.T.reshape((-1,9)) 

        #     extras["BEC"] = bec.tolist()

        #     # file = self.opts["save"]["BEC"]
        #     # if file is not None :
        #     #     with open(file, 'a') as f:
        #     #         line = "{:d} : {:s}\n".format(self.count,bec.flatten().tostring())
        #     #         f.write(line)
        #     #         # np.savetxt(f, bec.flatten(),delimiter=" ")

        return pot, force, vir, json.dumps(extras)
