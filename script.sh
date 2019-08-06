#!/bin/sh

set -- *.wiki
while [ $# -gt 0 ]
do
  pandoc -f creole -t markdown -s "${1}" -o "${1%.wiki}.md" &
  shift
  echo "pandoc -f creole -t markdown -s ${1} -o ${1%.wiki}.md"
  if [ $# -gt 0 ]
  then
    pandoc -f creole -t markdown -s "${1}" -o "${1%.wiki}.md" &
    shift
  fi
  wait
done
