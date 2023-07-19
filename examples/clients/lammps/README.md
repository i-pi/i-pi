 - client: Lammps (https://www.lammps.org)

   The command to inform lammps that it will run as client code is: 
       ```
       fix ID group-ID ipi address port [unix] [reset] 
       ```
 
   where
   ID, group-ID are documented in fix command
   ipi = style name of this fix command
   address = internet address (FQDN or IP), or UNIX socket name
   port = port number (ignored for UNIX sockets)
   keyword = unix or reset
       unix args = none = use a unix socket
       reset args = none = reset electrostatics at each call

   Check a more detailed documentation at https://docs.lammps.org/fix_ipi.html
