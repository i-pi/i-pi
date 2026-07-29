# Resonant electric field through an i-PI socket

This is the socket counterpart of
[`../mace-polar+E-resonant-client`](../mace-polar+E-resonant-client/README.md).
i-PI evaluates the resonant plane-wave field and sends its instantaneous value
to `extmace` for every force request.

Run it from this directory:

```bash
./run.sh
```

The script uses the shared resonant-client model assets, starts i-PI, and
connects `i-pi-py_driver` through the UNIX socket
`mace-polar-resonant-extra`. The driver is launched with
`requires_extra=true`, so each request follows:

```text
NEEDEXTRA → EXTRADATA(JSON with Efield) → READY → POSDATA → GETFORCE
```

The field amplitude is 0.1 V/Angstrom along z and its angular frequency is
0.719335584407 rad/fs. The final three columns of `properties.out` show the
time-dependent field in V/Angstrom.
