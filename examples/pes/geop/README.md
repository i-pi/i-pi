example with qTIP4P/f water model
---------------------------------

One runs i-PI as usual, followed by one or more instances of the driver:

```bash
i-pi input.xml &
i-pi-driver -u -h driver -m qtip4pf
```

Remember that the `<i-pi-root>/env.sh` file must be sourced before running
`i-pi` or `i-pi-driver`:

```bash
source <i-pi-root>/env.sh
```
