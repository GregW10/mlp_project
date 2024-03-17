import math
import sys


def b3s(sizeof_long_double: int, sizeof_flt: int, num_epochs: int) -> int:
    return 12 + 4*sizeof_long_double + 6*sizeof_flt*(1 + 3*num_epochs)


ld_size = 16

masses = {1.0: 1.65*(60/16.0), 100.0: 165.0*(60/16.0), 10_000.0: 16_500.0*(60/16.0),
          1_000_000.0: 1_650_000.0*(60/16.0),
          100_000_000.0: 165_000_000.0*(60/16.0)}

denom = 2**6

vel_scale_step = 1/denom

evolutions = 16
epochs = 1_000
iterations = 1_000_000
num_mass_means = evolutions # so that each velocity scaling factor has a large mass range
# starting_mass = 1.0 #variable
# mass_step = 1 #variable
mass_sd_scaling = 0.5
radius = 0
# timestep = 1/2**24 # use default
timestep_scaling = 1
timestep_exp = 0.25
# velocity_scaling = 0
# box_width = 2 # use default
# box_length = 2 # use default
# box_height = 2 # use default
# softening = 1/2**26 # use default
early_stopping = True
verbose = True
output_nsys = False
# threads = 60 # use default

soft = False


def main():
    print(f"Velocity scaling step = {vel_scale_step}")
    print(f"Masses = {masses}")
    num_sims = evolutions*len(masses)*(denom + 1)
    print(f"Total number of simulations = {num_sims}")
    single_fsize = b3s(ld_size, ld_size, epochs)
    print(f"Size of a single .3bod file = {single_fsize}")
    tot_fsize = num_sims*single_fsize
    print(f"Total .3bod file size of all simulations = {tot_fsize} bytes")
    counter = 0
    velocity_scaling = 0
    fname = "allgens.sh"
    f = open(fname, "w")
    if sys.platform.startswith("darwin"):
        f.write("#!/bin/zsh\n\n")
    elif sys.platform.startswith("linux"):
        f.write("#!/bin/bash\n\n")
    else:
        raise Exception("Error: invalid platform.\n")
    while velocity_scaling <= 1:
        print(f"vel scale = {velocity_scaling}")
        for starting_mass, mass_step in zip(masses.keys(), masses.values()):
            f.write(f"nohup gen -ve --evolutions {evolutions} --epochs {epochs} --iterations {iterations} "
                    f"--radius {radius} --num_mass_means {num_mass_means} "
                    f"--starting_mass {starting_mass} --mass_step {mass_step} --mass_sd_scaling {mass_sd_scaling} "
                    f"--timestep_scaling {timestep_scaling} --timestep_exp {timestep_exp} "
                    f"--velocity_scaling {velocity_scaling} {'--softening 0 ' if not soft else ''}"
                    f"-o vel_scale_{velocity_scaling}_sm_{starting_mass}_ms_{mass_step} > "
                    f"log_vel_scale_{velocity_scaling}_sm_{starting_mass}_ms_{mass_step} 2> "
                    f"errors_vel_scale_{velocity_scaling}_sm_{starting_mass}_ms_{mass_step}\n\n")
            f.write(f"echo 'vel_scale_{velocity_scaling}_sm_{starting_mass}_ms_{mass_step}' DONE\n\n")
            print(f"Starting mass = {starting_mass}, final mass = {starting_mass + (num_mass_means - 1)*mass_step}")
        counter += 1
        velocity_scaling += vel_scale_step
    f.close()
    print(f"Denom = {denom}")
    print(f"Counter = {counter}")


if __name__ == "__main__":
    main()
