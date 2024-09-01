# clkcomb

`clkcomb` is a program developed for the combination of multi-GNSS clock and
phase bias products. It is built on an earlier version of the
[PPPx](https://github.com/YuanxinPan/PPPx_bin) software package.

The initial development was completed as part of my MSc thesis in 2021. Since then,
the software has been further developed and applied to the IGS repro3 project at
Wuhan University. For a detailed explanation of the theories behind clock and
phase bias combination, please refer to my [MSc thesis](doc/Pan_MSc_Thesis_2021.pdf).

The following publications are based on this software (or an updated version of it):
- Geng J., Yan Z., Wen Q., Männel B., Masoumi S., Loyer S., Mayer-Gürr T., & Schaer S. (2024)
  Integrated satellite clock and code/phase bias combination in the third IGS
  reprocessing campaign. GPS Solut 28(3), 1-19.
  https://doi.org/10.1007/s10291-024-01693-9
- Geng J., Wen Q., Chen G., Dumitraschkewitz P., & Zhang Q. (2024)
  All-frequency IGS phase clock/bias product combination to improve PPP
  ambiguity resolution. J Geod 98(6), 48.
  https://doi.org/10.1007/s00190-024-01865-y


## Installation

### Linux & MacOS

Ensure you have the GNU make and a valid c++ compiler installed on your computer.
To install `clkcomb`, please clone the repository and compile the source code:

```
git clone git@github.com:YuanxinPan/clkcomb.git
cd clkcomb && make
```

The executables will be generated in the `bin/` directory.


### Windows

Open `src/vc++/clkcomb.sln` with Visual Studio (version 2013 or later).
Please note that the program has not yet been tested on Windows.


## Running the tests

An example is provided in the
[PPPx\_bin](https://github.com/YuanxinPan/PPPx_bin/tree/main/example/clkcomb)
repository to help you get started.


## Contributing

Please contribute via the **Pull requests** system.


## Author

- **Yuanxin Pan** - [YuanxinPan](https://github.com/YuanxinPan)


## License

This project is licensed under the BSD-3-Clause License - see the
[LICENSE](LICENSE) file for details.

