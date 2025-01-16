site=$1

valid_sites=("Grenoble" "Lille" "Luxembourg" "Lyon" "Nancy" "Nantes" "Rennes" "Sophia" "Strasbourg" "Toulouse")

if [ -z "$site" ]; then
    echo "Error: No site specified."
    echo "Valid sites are:"
    echo "${valid_sites[*]}"
    exit 1
fi

shopt -s nocasematch
if [[ ! " ${valid_sites[@]} " =~ " ${site} " ]]; then
    echo "Error: Invalid site specified."
    echo "Valid sites are:"
    echo "${valid_sites[*]}"
    exit 1
fi
shopt -u nocasematch

cd ~/Study/
tar --exclude='build'   \
    --exclude='.vscode' \
    --exclude='.git'    \
    -czvf reduction.tar.gz MPIReduction
scp reduction.tar.gz $site.g5k:
ssh $site.g5k "rm -rf MPIReduction"
ssh $site.g5k "tar -xzvf reduction.tar.gz"
ssh $site.g5k "cd MPIReduction && mkdir build && cd build && cmake .. && make"