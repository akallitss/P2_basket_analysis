# WELCOME TO THE ROMAPPING EXPERIENCE!
Mapping like you have never seen before...

## Introduction
This is a series of scrips that tries to automatise the process of creating mapping files from a Gerber file.
Please note that Gerber files are purely graphical representations and, as such, the correct setting for the following scripts might in principle be dependent on the "style" adopted by the designer (so careful cross-checking is highly advised).

To produce a mapping, given a comprehensive set of Gerber files, the following steps must be executed in sequence:

## Setp 1: Matching
The first script reads the Gerber files containing tracks and pins and produces a json file that links all pins to the **position** of the pad **connector** (NOT the pad itself!).<br>
To do this, one must first identify and provide the **pin file** (containing the graphic representation of the connectors) and any number of **track files** (containing the tracks leading from the pins to the pads).<br>
Subsequently, the user input will be required through a GUI to identigy the different **connectors** and the **pad area**.<br>
**NOTE:** The connectors require the right orientation to correctly order the pins! This is achievable through the <kbd>Rotate</kbd> and <kbd>Flip</kbd> buttons located on top of the GUI. Normally, a wrong orientation will result in missing pads in the following "pad selection" GUI but any valid configuration is always defined up to 180<sup>$\circ$</sup> connector rotations (recoverable later with the following scripts).

One also has the option to select multiple pad areas (say if you have pads with different pitches on the same board). In that case, if one area contains another, one should select first the innermost area since the code will then proceed to match the areas in order and remove any connections performed on the previous area from the latters.

A possible code snippet goes as follow:
```bash
python3 Matching.py --in_dir <gerber_folder> \
    --track_files <gerber_track_file_1> <gerber_track_file_2> ... \
    --pin_file <gerber_pin_file> \
    --n_connectors <number_of_connectors_to_map> \
    --n_pad_areas <number_of_pad_areas> \
    --checks #<-[OPTIONAL]: shows the connections right before making the pin selection
```
**[TO DO]:** Add zoom on the GUI, add -h option to the code

## Setp2: Discretization
Taking the previous output as input, it indexes the pad position
NOTE: Since the "pad position" is actually the position of the track-to-pad connection, this might not be the pad center and as such indexing might be non trivial (it is highly suggested to modify this file to index this according to an appropriate custom logic).
```bash
python3 Discretize.py --in_file <path_to_file><in_file>
```
**[TO DO]:** Add arguments for pad sizes and indexing logic, add-h option to the code

## Setp3(Optional): Transform
Taking the discretized mapping as input, it applies optional transformations to the pad cooedinates or pin indices and outputs a transformed file.
```bash
python3 Transform.py --in_file <pad_indexed_mapping> \
    --out_file <transformed_mapping> \
    --rotate_pads <n_90deg_rotations> \ #<-[OPTIONAL]: rotates the pad indices by 90<sup>$\circ$</sup> each time
    --flip_pads \ #<-[OPTIONAL]: flip pad x and y indices (consequently flipping the detector)
    --rotate_connectors #<-[OPTIONAL]: rotates the indexing of the connector
```

## Setp4: Format
Format the output four your preferred analysis software

### TBreco
```bash
python3 TBreco_format.py --in_file <path_to_file><pad_indexed_mapping> 
```
