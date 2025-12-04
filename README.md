 # GMT FEM transfer functions

A CLI interface to compute transfer functions between a set of inputs and outputs of the GMT finite element model (FEM).

The CLI is compiled for a specific FEM with the location of the model given by the environment variable `FEM_REPO` e.g.
```shell
export FEM_REPO=<path-to-FEM>
```

Contextual help is shown with
```shell
cargo r -r -- help
```
or for a specific sampling frequency command
```shell
cargo r -r -- help log-space
```
## Examples

Computing the transfer functions between M1 hardpoints and M1 & M2 rigid body motions, at the frequencies [1,5,10,20]Hz (the transfer functions are saved by default to the pickle file: `gmt_frequency_response.pkl`):
```shell
cargo r -r -- -i oss-harpoint-delta-f -o ossm1-lcl -o mcm2-lcl6-d set -v 1 -v 5 -v 10 -v 20
```

Computing the transfer functions between M1 hardpoints and the segment tip-tilt in the focal plane, logarithmically sampled with 1000 frequencies between 0.01Hz and 100Hz and saving the transfer functions to a Matlab mat file:
```shell
cargo r -r -- -i oss-harpoint-delta-f -o segment-tip-tilt -f m1-hp_segment-tt.mat log-space -l 0.01 -u 100 -n 1000
```

# Installation

Instead of running from the crate location, a executable binary can be compiled locally with:

```shell
FEM_REPO=<path-to-FEM> cargo install --git https://github.com/rconan/gmt-fem-transfer-functions.git --bin gmt-fem-transfer-functions
```

The binary is installed locally and the CLI app is called simply by invoking `gmt-fem-transfer-functions` (but still with `FEM_REPO` set to the path the GMT FEM) e.g.
```shell
FEM_REPO=<path-to-FEM> gmt-fem-transfer-functions -i oss-harpoint-delta-f -o segment-tip-tilt -f m1-hp_segment-tt.mat log-space -l 0.01 -u 100 -n 1000
```

