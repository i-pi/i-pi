The example folders is structured such that each sub-folder is focused on a different aspect of using i-PI:

- **clients**  :    Contains examples that are code-specific, highlighting how the driver code should be set up
                    (client-specific syntax and tags) to run it properly with i-PI

- **features** :    Examples of different functionalities implemented in i-PI. 
                    All examples can be run locally with the drivers provided with the code.

- **hpc_scripts** :    Examples of submission scripts on HPC platforms

- **temp**     :    Temporary folder with historic examples that have not yet been adapted
                    to the current folder structure

- **init_files**: repository of input files shared by many examples

We keep this folder updated as much as we can, and try to run automated tests on these inputs, but in some cases, e.g. when using external clients, we cannot run tests. 
Please report a bug if you find something that is not working. 

If you want to add an example for a new feature, please respect the folder structure above. In addition, please:
- Keep the size of the example to the minimum. If you need large files to run the client (e.g. a pseudopotential, or the parameters of a ML potential) give instructions on how to download the file from a repository. If the file is shared by many examples, put it in `init_files` and make a symlink
- Include a README that briefly explains what the example is about, and what a user should be aware of when running it (what it does, what special setup is needed, etc.)

Notes on file naming conventions and formats:
  - Please use lowercase names rather than Camel_Case
  - Please use use underscores rather than hyphen 
  - Please use driver whenever possible
