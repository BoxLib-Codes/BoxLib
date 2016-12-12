#!/usr/bin/env bash

OLD_SRC_PATH="Tools\/C_mk"
NEW_SRC_PATH="Tools\/GNUMake"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -type d \( -name .git -o -path './Tools/Migration' \) -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

echo "Remove 'FCOMP = ' line from GNUmakefile"
find . -type d \( -name .git -o -path './Tools/Migration' \) -prune -o -type f -name "GNUmakefile" -exec sed -i '/FCOMP\s*=/d' {} +
