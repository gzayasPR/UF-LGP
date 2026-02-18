#!/bin/bash

proj_dir="$(pwd)"
my_bin="$(pwd)/bin"
my_metadata="$(pwd)/metadata"
my_results="$(pwd)/results"
my_docs="$(pwd)/docs"
my_bash="$(pwd)/bash_out"
mkdir -p $my_bin
mkdir -p $my_metadata
mkdir -p $my_results
mkdir -p $my_docs
mkdir -p $my_bash
env_file="$my_bin/project_env.sh"

cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="$(pwd)"
export my_bin="$(pwd)/bin"
export my_metadata="$(pwd)/metadata"
export my_results="$(pwd)/results"
export my_docs="$(pwd)/docs"
export my_bash="$(pwd)/bash_out"
EOL

# Make the file executable
chmod +x "$env_file"


cd $my_bin

