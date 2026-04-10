#!/bin/bash

GRO=$1
TPR=$2

# replace extension with _mol.<ext>
MOL_NAME=$(echo "$GRO" | sed 's/\.[^.]*$/_mol&/')
# replace extension with _center.<ext>
CENTER_NAME=$(echo "$GRO" | sed 's/\.[^.]*$/_center&/')

# Make molecules whole and centre on Protein; output System
echo -e "1\n0" | gmx trjconv -f "$GRO" -s "$TPR" -pbc mol -center -o "$MOL_NAME"

# Re-centre in compact box; centre on Protein, output System
echo -e "1\n0" | gmx trjconv -f "$MOL_NAME" -s "$TPR" -center -pbc mol -boxcenter zero -ur compact -o "$CENTER_NAME"
