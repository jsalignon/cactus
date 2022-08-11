#!/bin/bash

# cd /home/jersal/workspace/cactus/software
# docs/menu/add_menu.sh

menu_file=docs/menu/menu.md
# doc_files=$(ls docs/*.md)
# doc_files=$(find docs -name "*.md")
# readarray -d '' doc_files < <(find docs -name "*.md" -not -path "docs/menu/*" -print0)
# echo ${doc_files[@]} | grep READ
IFS=$'\n'
doc_files=($(find docs -name "*.md" -not -path "docs/menu/*"))
unset IFS
doc_files+=("README.md")
# echo ${doc_files[@]} | grep READ
# doc_file=${doc_files[11]}

for doc_file in ${doc_files[@]}
do
  awk -i inplace 'x==1 {print} /END_OF_MENU/ {x=1}' $doc_file
  echo "$(cat $menu_file <(echo) <(echo) $doc_file)" > $doc_file
done

