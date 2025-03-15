#!/bin/bash
#
# create_sim.sh
#
# Purpose:
#   1) Prompt user for:
#      - Simulation name
#      - A list of salts (user-input)
#      - Number of data files N per salt
#      - Data folders & percentages (must sum to 100)
#   2) For each salt in the user-defined list:
#      - For each data folder, randomly pick (percentage*N/100) from data/<folder>/<salt>/new_data_*
#      - Copy them to SIMNAME/<salt>, renamed as new_data_1..new_data_k
#      - Copy snippet files, do sed replacements in create_slab.in & run_slab_mixture.in:
#         "s/EECC/$eecc/", "s/SSCC/$sscc/", "s/XXXX/$rndm/"
#      - eecc, sscc read from salt_map_E.py, salt_map_A.py
#
# Run from inside slabs/:  ./create_sim.sh
#

###############################################################################
# 1) Ask user for simulation name, number of salts, etc.
###############################################################################

read -p "Enter simulation name: " SIMNAME
mkdir -p "$SIMNAME"
echo "Created folder: $SIMNAME"

# Ask how many salts, then read each salt name into an array
read -p "How many salts do you have? " NSALTS
declare -a SALTS
for (( i=1; i<=NSALTS; i++ )); do
  read -p "Enter salt #$i (e.g. 0.073): " S
  SALTS+=( "$S" )
done

read -p "How many copies of protein (N) per salt? " N

read -p "How many types to mix? " NFOLDERS

declare -a FOLDERS
declare -a PERCENTAGES
sum_perc=0

for (( i=1; i<=NFOLDERS; i++ )); do
  echo
  echo "=== Data folder $i of $NFOLDERS ==="
  read -p "Folder name (e.g. CENPB-172, Control): " fname
  read -p "What percentage of $N should come from '$fname'? " perc

  FOLDERS+=( "$fname" )
  PERCENTAGES+=( "$perc" )
  (( sum_perc += perc ))
done

if (( sum_perc != 100 )); then
  echo "ERROR: Percentages must sum to 100. Got $sum_perc."
  exit 1
fi

###############################################################################
# 2) We'll keep each folder+salt's new_data_* in data_arrays["folder:salt"].
###############################################################################

declare -A data_arrays

init_data_array() {
  local folder="$1"
  local salt="$2"
  local path="data/$folder/$salt"

  if [[ -d "$path" ]]; then
    local all_files
    all_files=$(find "$path" -type f -name "new_data_*" 2>/dev/null)
    data_arrays["$folder:$salt"]="$all_files"
  else
    # folder/salt subdir not found => no files
    data_arrays["$folder:$salt"]=""
  fi
}

draw_random_files() {
  local folder="$1"
  local salt="$2"
  local count="$3"

  local key="$folder:$salt"
  local list="${data_arrays["$key"]}"

  # Convert multiline string to array
  readarray -t array <<< "$list"
  local total="${#array[@]}"

  if (( total < count )); then
    echo "ERROR: $folder/$salt only has $total files, need $count."
    exit 1
  fi

  # Shuffle
  mapfile -t shuffled < <(printf "%s\n" "${array[@]}" | shuf)

  # First 'count' are selected
  local selected=("${shuffled[@]:0:count}")
  local remain=("${shuffled[@]:count}")

  # Update the array with leftover
  data_arrays["$key"]=$(printf "%s\n" "${remain[@]}")

  # Output the selected
  printf "%s\n" "${selected[@]}"
}

###############################################################################
# 3) For each salt, pick the fraction from each folder, rename them, copy snippet
###############################################################################

bash create_lammps_code.sh

# The snippet's "files to copy"
files_to_copy=(
  "data_files.in" 
  "chromatin.so"
  "NAFlex_params.txt"
  "create_slab.in"
  "run_slab.sub"
  "run_slab_mixture.in"
  "run_slab.sh"
  "create_slab.sub"
  "DNA_sequence.txt"
)

files_to_copy_parent=(
  "chromatin.so"
  "NAFlex_params.txt"
  "density.py"
  "submit_all_create.sh"
  "phase_diagram.py"
  "run_slab.sh"
)


# We'll do "s/XXXX/$rndm/" as requested:
rndm=$RANDOM
XXX=192673  # As per your snippet, though unused directly in sed

echo
echo "===== Distributing data across salts, applying snippet logic ====="

for salt in "${SALTS[@]}"; do
  echo
  echo "Processing salt: $salt"

  # Pre-initialize data_arrays for each (folder, salt)
  for (( i=0; i<NFOLDERS; i++ )); do
    folder="${FOLDERS[$i]}"
    key="$folder:$salt"
    if [[ -z "${data_arrays["$key"]}" ]]; then
      init_data_array "$folder" "$salt"
    fi
  done

  # We'll gather all picks for this salt in an array
  chosen=()
  total_chosen=0


  # Gather all available files from all folders
  all_files=()
  for folder in "${FOLDERS[@]}"; do
    key="$folder:$salt"
    if [[ -z "${data_arrays["$key"]}" ]]; then
      init_data_array "$folder" "$salt"
    fi
    # Append all files from this folder to a global list
    mapfile -t files <<< "${data_arrays["$key"]}"
    all_files+=( "${files[@]}" )
  done

  # Shuffle the combined file list
  mapfile -t shuffled < <(printf "%s\n" "${all_files[@]}" | shuf)

  # Pick the first N files from the shuffled list
  chosen=( "${shuffled[@]:0:N}" )
  total_chosen="${#chosen[@]}"

  echo "  -> Picked $total_chosen mixed files for salt '$salt'."

  #for (( i=0; i<NFOLDERS; i++ )); do
  #  folder="${FOLDERS[$i]}"
  #  perc="${PERCENTAGES[$i]}"

  #  count=$(( (perc * N) / 100 ))
  #  if (( count > 0 )); then
  #    # Randomly pick 'count' from data/<folder>/<salt>
  #    mapfile -t picks < <(draw_random_files "$folder" "$salt" "$count")
  #    chosen+=( "${picks[@]}" )
  #    (( total_chosen += count ))
  #  fi
  #done

  #echo "  -> Picked $total_chosen files for salt '$salt' (ideal was $N)."

  # Create the salt directory under SIMNAME
  dir="$SIMNAME/$salt"
  mkdir -p "$dir"
  mkdir "$dir/data"

  # Rename each chosen file new_data_1.. up to total_chosen
  for (( i=0; i<total_chosen; i++ )); do
    old_path="${chosen[$i]}"
    new_id=$(( i + 1 ))
    new_file="new_data_$new_id"
    cp "$old_path" "$dir/data/$new_file.txt"
  done

  # Copy snippet files
  for f in "${files_to_copy[@]}"; do
    cp "$f" "$dir/"
  done

  for f in "${files_to_copy_parent[@]}"; do
    cp "$f" "$SIMNAME/"
  done
  
  # Get eecc, sscc from Python
  eecc=$(python3 salt_map_E.py "$salt")
  sscc=$(python3 salt_map_A.py "$salt")
  echo "  -> eecc=$eecc, sscc=$sscc"

  # Perform sed replacements in create_slab.in & run_slab_mixture.in
  sed -i "s/EECC/$eecc/" "$dir/create_slab.in"
  sed -i "s/EECC/$eecc/" "$dir/run_slab_mixture.in"

  sed -i "s/SSCC/$sscc/" "$dir/create_slab.in"
  sed -i "s/SSCC/$sscc/" "$dir/run_slab_mixture.in"

  sed -i "s/XXXX/$rndm/" "$dir/run_slab_mixture.in"
done

echo
echo "Done! Each salt in '$SIMNAME' has up to N data files, plus snippet files with sed replacements."
