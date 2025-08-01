"""Small functions/classes providing access to driver PES to be called from driver.py"""

import pkgutil
import importlib
import traceback

__all__ = []

# Dictionary to store driver name to class mapping
__drivers__ = {}

# Iterate through all modules in the current package folder
for loader, module_name, is_pkg in pkgutil.iter_modules(__path__):
    # Import the module
    try:
        module = importlib.import_module("." + module_name, __package__)
    except Exception:
        print(f"!! Could not import PES module {module_name} !!")
        traceback.print_exc()

    # Get the driver class and name from the module
    driver_class = getattr(module, "__DRIVER_CLASS__", None)
    driver_name = getattr(module, "__DRIVER_NAME__", None)

    # If both class and name are defined, update __all__ and __drivers__
    if driver_class and driver_name:
        __all__.append(driver_class)
        __drivers__[driver_name] = getattr(module, driver_class)
        globals()[driver_class] = getattr(module, driver_class)  # add class to globals
    else:
        if driver_class != "driver_tools":
            raise ImportError(
                f"PES module `{module_name}` does not define __DRIVER_CLASS__ and __DRIVER_NAME__"
            )

__all__.append("__drivers__")


def gpu_oversubscription():
    """
    Distributes GPUs among SLURM tasks by oversubscribing them when
    there are more tasks than available GPUs on a node.

    Determines the GPU assignment based on SLURM_LOCALID and the total
    number of GPUs and tasks, sets CUDA_VISIBLE_DEVICES accordingly,
    and logs output to a file named 'task.localid={ID}.log'.
    """

    import os, sys

    # Get the task's local ID on this node (used for GPU assignment and logging)
    local_id = os.environ.get("SLURM_LOCALID", "unknown")

    # Redirect stdout and stderr to a task-specific log file
    log_file = f"task.localid={local_id}.log"
    sys.stdout = open(log_file, "w")
    sys.stderr = sys.stdout  # Also capture warnings and errors

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
