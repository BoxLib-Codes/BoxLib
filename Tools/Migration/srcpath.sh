#!/usr/bin/env bash

OLD_SRC_PATH="C_AMRLib"
NEW_SRC_PATH="Amr"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -path './.git' -prune -o -path './Tools/Migration' -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_AmrCoreLib"
NEW_SRC_PATH="AmrCore"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -path './.git' -prune -o -path './Tools/Migration' -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_BaseLib"
NEW_SRC_PATH="Base"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -path './.git' -prune -o -path './Tools/Migration' -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_BoundaryLib"
NEW_SRC_PATH="Boundary"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -path './.git' -prune -o -path './Tools/Migration' -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="C_ParticleLib"
NEW_SRC_PATH="Particle"
echo ${OLD_SRC_PATH}" --> "${NEW_SRC_PATH}
find . -path './.git' -prune -o -path './Tools/Migration' -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +


