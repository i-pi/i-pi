# Static electric field through an i-PI socket

This is the socket counterpart of
[`../mace-polar+E-client`](../mace-polar+E-client/README.md). i-PI evaluates a
static 1 V/Angstrom field along z and sends it to `extmace` before each force
request.

Run it from this directory:

```bash
./run.sh
```

The script downloads the shared MACE-POLAR checkpoint if necessary, launches
i-PI, and then launches `i-pi-py_driver`. Its `requires_extra=true` parameter
makes the driver advertise `NEEDEXTRA`; i-PI responds with an `EXTRADATA`
message containing UTF-8 JSON such as:

```json
{"time": 0.0, "where": "client", "Efield": [0.0, 0.0, 0.0194469]}
```

`extmace` parses that JSON through `store_extra()` and applies `Efield` to the
MACE calculation. The field values in the example XML are converted to atomic
units before transmission. `properties.out` prints the same field in V/Angstrom.
