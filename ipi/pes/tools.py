import os
import time
from datetime import datetime


# --------------------------------------- #
def gpu_oversubscription():
    """
    Assigns GPUs to SLURM tasks by oversubscribing when tasks exceed GPUs.

    Uses SLURM_LOCALID and SLURM_NTASKS_PER_NODE to map tasks onto GPUs,
    sets CUDA_VISIBLE_DEVICES accordingly, and prints assignment info.
    """

    # enable this in your submission script by using
    # export IPI_DRIVER_GPU=true
    # which will turn on the GPU oversubscription
    if os.environ.get("IPI_DRIVER_GPU", False):

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


# --------------------------------------- #
class Timer:
    _stack = []
    _records = []

    def __init__(self, name, log=True):
        self.name = name
        self.log = log
        self.level = len(Timer._stack)

    def __enter__(self):
        self.start = time.time()
        Timer._stack.append(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        self.elapsed = self.end - self.start
        Timer._stack.pop()
        if self.log:
            Timer._records.append((self.start, self.level, self.name, self.elapsed))

    @staticmethod
    def report():
        for t, level, name, elapsed in sorted(Timer._records):
            timestamp = datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S")
            indent = "    " * level
            print(f"{timestamp}{indent} {name}: {elapsed:.4f} s")

    @staticmethod
    def reset():
        Timer._records.clear()
