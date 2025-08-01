"""An interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

import json

from .ase import ASEDriver

MACECalculator = None

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"


class MACE_driver(ASEDriver):
    """
    Driver for the MACE MLIPs.
    The driver requires specification of a torch model,
    and a template file that describes the chemical makeup
    of the structure.

    Command-line:
    i-pi-py_driver -m mace -u -a address -o template=structure.xyz,model=mace_mp.model

    Parameters:
    :param template: string, filename of an ASE-readable structure file
        to initialize atomic number and types
    :param model: string, filename of the MACE model
    """

    def __init__(
        self,
        template,
        model,
        device="cpu",
        requires_extra: bool = False,
        gos: bool = False,  # GPU OverSubscription
        mace_kwargs=None,
        *args,
        **kwargs,
    ):
        # warning(
        #     "THIS PES HAS NOT BEEN TESTED FOLLOWING CONVERSION TO THE NEW PES API.",
        #     verbosity.low,
        # )
        global MACECalculator

        try:
            from mace.calculators import MACECalculator
        except:
            raise ImportError("Couldn't load mace bindings")

        try:
            from ase.outputs import _defineprop, all_outputs

            # avoid duplicate
            # it complains with a committee of ffdirect MACE models
            if "node_energy" not in all_outputs:
                _defineprop("node_energy", dtype=float, shape=("natoms",))
        except ImportError:
            raise ValueError("Could not find or import the ASE module")

        self.model = model
        self.device = device
        self.mace_kwargs = {}
        if mace_kwargs is not None:
            with open(mace_kwargs, "r") as f:
                self.mace_kwargs = json.load(f)

        super().__init__(template, *args, **kwargs)
        self.requires_extra = requires_extra

        if gos:
            gpu_oversubscription()

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_parameters()

        self.ase_calculator = MACECalculator(
            model_paths=self.model, device=self.device, **self.mace_kwargs
        )


def gpu_oversubscription():
    """
    Assigns GPUs to SLURM tasks by oversubscribing when tasks exceed GPUs.

    Uses SLURM_LOCALID and SLURM_NTASKS_PER_NODE to map tasks onto GPUs,
    sets CUDA_VISIBLE_DEVICES accordingly, and prints assignment info.
    """

    import os

    # Get the task's local ID on this node (used for GPU assignment and logging)
    local_id = os.environ.get("SLURM_LOCALID", "unknown")

    print("\n[GPU Oversubscription Info]")

    # Get total number of GPUs available on this node
    num_gpus = len([gpu for gpu in os.popen("nvidia-smi -L").readlines()])
    print(f"  - Number of GPUs available: {num_gpus}")

    # Parse the current task's local ID (default to 0 if not set)
    local_id = int(os.environ.get("SLURM_LOCALID", 0))
    print(f"  - SLURM_LOCALID: {local_id}")

    # Total number of tasks assigned to this node
    num_tasks = int(os.environ.get("SLURM_NTASKS_PER_NODE", 1))
    print(f"  - SLURM_NTASKS_PER_NODE: {num_tasks}")

    # Compute how many tasks will share each GPU
    tasks_per_gpu = max(1, num_tasks // max(num_gpus, 1))
    print(f"  - Estimated tasks per GPU: {tasks_per_gpu}")

    # Assign a GPU ID to this task based on its local ID
    gpu_id = local_id // tasks_per_gpu
    gpu_id = min(gpu_id, num_gpus - 1)  # Clamp to last available GPU
    print(f"  - Assigned GPU ID: {gpu_id}")

    # Restrict visibility to the assigned GPU
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    print(f"  - CUDA_VISIBLE_DEVICES set to: {os.environ['CUDA_VISIBLE_DEVICES']}")
