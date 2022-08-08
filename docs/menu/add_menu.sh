#!/bin/bash

# cd /home/jersal/workspace/cactus/software
# docs/menu/add_menu.sh

menu_file=docs/menu/menu.md
doc_files=$(ls docs/*.md)
doc_files+=("README.md")

for doc_file in ${doc_files[@]}
do
  sed -i 1,5d $doc_file
  echo "$(cat $menu_file <(echo) $doc_file)" > $doc_file
done

