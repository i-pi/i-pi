These examples show how connect i-PI to client code using  ASE as a middle man.
There are two possible ways of create that connection.

a) i-PI (server) <----> ASE(client) <---> Force code
b) i-PI (server) <----> ASE(client) <----> ASE(server) <---> Force code (client)

While the first option restarts the code for any new calculation, the second one maintains it alive but
might be a little bit trickier to setup in a  cluster.
