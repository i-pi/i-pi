Have a look at the ```README.md``` in ```../efield-water/```

**Pay attention** that ```bec.txt``` **must** be modified accordingly to the number of ```beads```!
If you have ```beads=N``` and you want to have the same BEC tensors for all the beads, you just need to repeat the same array ```N``` times in ```bec.txt```.
This means that in ```bec.txt``` you should provide an array of shape ```(nbeads x 3 x Natoms, 3)```.
You can provide the same array (flattened) directly in the ```input.xml``` file.
