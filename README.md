This is a plugin for [PSI4](https://github.com/psi4/psi4public) that writes the natural orbitals
to a molden file for visualization. It's been tested on 9cb03c3ff (It will NOT work in 4.0b5).

You will need a C++11 compiler (GCC 4.8 or Clang 3.4 will do). First generate a new makefile
using: psi4 --new-plugin-makefile and compile the plugin. There is an example `input.dat` to 
show how you can use the plugin.

The plugin doesn't calculate the naturals itself, you need to tell PSI4 to do that (using
the detci module for example). This plugin will read the calculate naturals and write them
out as a molden input file.

I've used this plugin for my own purposes. It works for me but your mileage may vary. 
I've tried to put a couple of lines of comments in the code to explain how to works.
